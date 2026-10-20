"""Deterministic tool layer: config, rate-limited clients, format scope, assay classification.

Nothing in this module calls a model. Every function here is a tool the orchestrator can call,
and every threshold is config, not a constant buried in a prompt.
"""
from __future__ import annotations

import csv, io, json, os, re, time, threading
from dataclasses import dataclass, field, asdict
from collections import Counter, defaultdict

import requests

# --------------------------------------------------------------------- config ---


@dataclass
class Config:
    # scope
    searchable_formats: tuple = ('.mzml', '.mzxml')
    convertible_formats: tuple = ('.raw', '.d')
    # packaging only. NOT .scan/.scans: those are index sidecars, i.e. a separate file that
    # describes an acquisition, never the acquisition itself.
    wrappers: tuple = ('.zip', '.gz', '.bz2', '.tar')
    sidecar_suffixes: tuple = ('.scan', '.scans', '.index', '.idx')
    # when the same acquisition is published twice (peak/ and ccms_peak/), prefer this order
    collection_preference: tuple = ('ccms_peak', 'peak', 'raw')
    # thresholds
    grammar_coverage: float = 0.95
    join_match_rate: float = 0.80
    join_injectivity: float = 1.05
    header_sample_per_subtree: int = 2
    # pacing (seconds between calls, per host)
    pace: dict = field(default_factory=lambda: {
        'massive.ucsd.edu': 12.0,              # ProteoSAFe 429s aggressively
        'metabolomics-usi.gnps2.org': 7.0,     # documented 10 req/min
        'datasetcache.gnps2.org': 1.0,
        'redu.gnps2.org': 1.0,
    })
    retries: int = 5
    cache_dir: str = '.agent_cache'

    def save(self, path):
        json.dump(asdict(self), open(path, 'w'), indent=1, default=list)


CONVERTERS = {
    ('.raw', 'thermo'): 'MSConvert',
    ('.raw', 'waters'): 'MSConvert (Symphony for MSe/SONAR/HDMSe/HD-DDA)',
    ('.raw', ''): 'MSConvert',
    ('.d', 'agilent'): 'MSConvert (MassHunter scan-number export has known issues)',
    ('.d', 'bruker'): 'CompassXport, or MSConvert when timsTOF (analysis.tdf present)',
    ('.d', ''): 'MSConvert or CompassXport depending on vendor',
}

PROTEOMICS_ONLY = ('.mzid', '.pepxml', '.fasta', '.fa', '.msf', '.dat', '.prot.xml', '.group')

# -------------------------------------------------------------------- http ------


class Http:
    """One rate-limited, retrying HTTP client. Treats a 500 page containing 429 as a 429."""

    def __init__(self, cfg: Config):
        self.cfg = cfg
        self._last = defaultdict(float)
        self._lock = threading.Lock()
        os.makedirs(cfg.cache_dir, exist_ok=True)

    def _wait(self, host):
        gap = self.cfg.pace.get(host, 0.5)
        with self._lock:
            delta = time.time() - self._last[host]
            if delta < gap:
                time.sleep(gap - delta)
            self._last[host] = time.time()

    def get(self, url, headers=None, timeout=240, want='text'):
        host = re.sub(r'^https?://([^/]+).*$', r'\1', url)
        last = None
        for attempt in range(self.cfg.retries):
            self._wait(host)
            try:
                r = requests.get(url, headers=headers or {}, timeout=timeout)
            except Exception as e:                                   # network flake
                last = ('exception', str(e)[:200])
                time.sleep(5 * (attempt + 1))
                continue
            body = r.text if want != 'bytes' else ''
            throttled = r.status_code == 429 or (r.status_code >= 500 and '429' in body[:4000])
            if throttled:
                time.sleep(45)
                last = ('throttled', r.status_code)
                continue
            if r.status_code >= 500:
                last = ('server', r.status_code)
                time.sleep(5 * (attempt + 1))
                continue
            return r
        raise RuntimeError('GET failed after %d attempts: %s (%s)' % (self.cfg.retries, url, last))

    def cached(self, url, name, **kw):
        """Fetch once per run; the cache is what makes a re-run cheap and auditable."""
        path = os.path.join(self.cfg.cache_dir, name)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            return open(path, encoding='utf-8', errors='replace').read()
        txt = self.get(url, **kw).text
        open(path, 'w', encoding='utf-8', newline='').write(txt)
        return txt


# ------------------------------------------------------------------ clients ----


class MassIVE:
    PROXI = 'https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/%s'
    DL = ('https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile'
          '?forceDownload=true&file=f.%s/%s')

    def __init__(self, http: Http):
        self.http = http

    def proxi(self, acc):
        return json.loads(self.http.cached(self.PROXI % acc, 'proxi_%s.json' % acc))

    def fetch_bytes(self, acc, path, nbytes=70000):
        r = self.http.get(self.DL % (acc, path), headers={'Range': 'bytes=0-%d' % nbytes},
                          timeout=300, want='bytes')
        if r.status_code == 200 and int(r.headers.get('Content-Length', 0)) > 5_000_000:
            raise RuntimeError('server ignored Range; refusing to pull %s bytes'
                               % r.headers.get('Content-Length'))
        return r.content[:nbytes]


class DatasetCache:
    BASE = 'https://datasetcache.gnps2.org/datasette/database'

    def __init__(self, http: Http):
        self.http = http

    def _table(self, table, acc):
        url = '%s/%s.csv?_size=max&dataset__exact=%s' % (self.BASE, table, acc)
        txt = self.http.cached(url, '%s_%s.csv' % (table, acc))
        return list(csv.DictReader(io.StringIO(txt)))

    def filenames(self, acc):
        return self._table('filename', acc)

    def uniquemri(self, acc):
        return self._table('uniquemri', acc)

    def sql(self, query, timeout=280):
        import urllib.parse
        url = '%s.json?_shape=array&sql=%s' % (self.BASE, urllib.parse.quote(query))
        r = self.http.get(url, timeout=timeout)
        return r.json() if r.status_code == 200 else []


class ReDU:
    """ReDU itself is the authority on its own schema - never hardcode the column list."""

    DUMP = 'https://redu.gnps2.org/dump'
    ATTRS = 'https://redu.gnps2.org/attributes'

    def __init__(self, http: Http):
        self.http = http

    def schema(self):
        """The 57-column header, read with a Range request instead of a ~930 MB download."""
        path = os.path.join(self.http.cfg.cache_dir, 'redu_schema.json')
        if os.path.exists(path):
            return json.load(open(path))
        r = self.http.get(self.DUMP, headers={'Range': 'bytes=0-8192'}, timeout=180)
        cols = r.text.split('\n')[0].rstrip('\r').split('\t')
        if len(cols) < 10 or cols[0] != 'ATTRIBUTE_DatasetAccession':
            raise RuntimeError('unexpected ReDU dump header: %r' % cols[:5])
        json.dump(cols, open(path, 'w'), indent=1)
        return cols

    def cv_columns(self):
        txt = self.http.cached(self.ATTRS, 'redu_attributes.json')
        return [a['attributename'] for a in json.loads(txt)]


class OLS4:
    """EBI Ontology Lookup Service v4 - https://www.ebi.ac.uk/ols4/api-docs

    Ontology identifiers are a lookup, never a recall task: the agent asks OLS4 for the term and
    its lineage. The model's only role here is picking between candidates OLS4 returns.
    """

    BASE = 'https://www.ebi.ac.uk/ols4/api'
    OBO = 'http://purl.obolibrary.org/obo/%s'

    def __init__(self, http: Http):
        self.http = http
        self._rank = {}

    def search(self, q, ontology=None, exact=False, rows=8):
        import urllib.parse
        p = {'q': q, 'rows': rows}
        if ontology:
            p['ontology'] = ontology
        if exact:
            p['exact'] = 'true'
        url = self.BASE + '/search?' + urllib.parse.urlencode(p)
        r = self.http.get(url, timeout=120)
        docs = r.json().get('response', {}).get('docs', []) if r.status_code == 200 else []
        return [{'obo_id': d.get('obo_id'), 'label': d.get('label'),
                 'ontology': d.get('ontology_name'), 'iri': d.get('iri'),
                 'description': (d.get('description') or [''])[0][:160]} for d in docs if d.get('obo_id')]

    def resolve(self, label, ontology):
        """Exact match first, then a loose search. Returns (chosen_or_None, candidates)."""
        hits = self.search(label, ontology=ontology, exact=True)
        if len(hits) == 1:
            return hits[0], hits
        loose = self.search(label, ontology=ontology, exact=False)
        exact_label = [h for h in loose if (h['label'] or '').lower() == label.lower()]
        if len(exact_label) == 1:
            return exact_label[0], loose
        return None, (hits or loose)

    def term(self, ontology, obo_id):
        iri = self.OBO % obo_id.replace(':', '_')
        import urllib.parse
        url = '%s/ontologies/%s/terms?%s' % (self.BASE, ontology, urllib.parse.urlencode({'iri': iri}))
        r = self.http.get(url, timeout=120)
        if r.status_code != 200:
            return None
        terms = r.json().get('_embedded', {}).get('terms', [])
        return terms[0] if terms else None

    def _rank_label(self, rank_iri):
        if rank_iri not in self._rank:
            short = rank_iri.rsplit('/', 1)[-1].replace('_', ':')
            t = self.term('taxrank', short)
            self._rank[rank_iri] = (t or {}).get('label')
        return self._rank[rank_iri]

    def lineage(self, taxid):
        """{rank label -> taxon label} for a taxon and its ancestors, ranks resolved via TAXRANK."""
        import urllib.parse
        obo = 'NCBITaxon:%s' % taxid
        own = self.term('ncbitaxon', obo)
        if not own:
            return {}, None
        iri = self.OBO % ('NCBITaxon_%s' % taxid)
        dbl = urllib.parse.quote(urllib.parse.quote(iri, safe=''), safe='')
        url = '%s/ontologies/ncbitaxon/terms/%s/hierarchicalAncestors?size=100' % (self.BASE, dbl)
        r = self.http.get(url, timeout=180)
        anc = r.json().get('_embedded', {}).get('terms', []) if r.status_code == 200 else []
        out = {}
        for t in anc + [own]:
            rank_iris = (t.get('annotation') or {}).get('has_rank') or []
            if not rank_iris:
                continue
            label = self._rank_label(rank_iris[0])
            if label and t.get('label'):
                out[label] = t['label']
        own_rank_iris = (own.get('annotation') or {}).get('has_rank') or []
        own_rank = self._rank_label(own_rank_iris[0]) if own_rank_iris else None
        return out, own_rank


class UsiResolver:
    URL = 'https://metabolomics-usi.gnps2.org/json/?usi1=%s'

    def __init__(self, http: Http):
        self.http = http

    def spectrum(self, usi, scan):
        import urllib.parse
        u = urllib.parse.quote('%s:scan:%d' % (usi, scan), safe=':/')
        r = self.http.get(self.URL % u, timeout=180)
        if r.status_code == 200 and r.text.startswith('{'):
            try:
                return r.json()
            except ValueError:
                return None
        return None

    def probe(self, usi, scans=(2, 50, 400, 1200)):
        """Returns (resolves, ms2_found, details). precursor_mz > 0 means an MS2 scan."""
        out = []
        for s in scans:
            d = self.spectrum(usi, s)
            if d:
                out.append({'scan': s, 'n_peaks': d.get('n_peaks'),
                            'precursor_mz': d.get('precursor_mz'),
                            'precursor_charge': d.get('precursor_charge')})
        return bool(out), any((o.get('precursor_mz') or 0) > 0 for o in out), out


# -------------------------------------------------------------------- scope ----


def strip_wrappers(name, cfg: Config):
    """Remove packaging suffixes, PRESERVING case.

    Case matters downstream: filename grammars are case-sensitive, so lowercasing here turns
    QC_mix_1 into qc_mix_1 and silently breaks every proposed regex.
    """
    b = name.replace('\\', '/').rstrip('/').split('/')[-1]
    changed = True
    while changed:
        changed = False
        for w in cfg.wrappers:
            if b.lower().endswith(w):
                b, changed = b[:-len(w)], True
    return b


def extension(name, cfg: Config):
    b = strip_wrappers(name, cfg).lower()
    if b.endswith('.mzdata.xml'):
        return '.mzdata.xml'
    m = re.search(r'(\.[a-z0-9]+)$', b)
    return m.group(1) if m else ''


def is_sidecar(name, cfg: Config):
    b = name.replace('\\', '/').rstrip('/').split('/')[-1].lower()
    return b.endswith(tuple(cfg.sidecar_suffixes))


def acquisitions(paths, cfg: Config):
    """Collapse vendor directories, drop sidecars, dedupe re-publications.

    Returns (acqs, census). An acquisition is a thing that was injected once, so:
      - every file inside a .d/.raw directory counts as that one directory,
      - a .mzML.scans / .wiff.scan index file is NOT an acquisition,
      - the same stem published under peak/ and ccms_peak/ is one acquisition, not two.
    """
    vendor_dirs, flat = {}, []
    for p in paths:
        m = re.match(r'^(.*?[^/]+\.(?:d|raw))/', p, re.I)
        if m:
            vendor_dirs[m.group(1)] = vendor_dirs.get(m.group(1), 0) + 1
        else:
            flat.append(p)

    census, candidates = Counter(), []
    for d in vendor_dirs:
        ext = extension(d, cfg)
        census[ext] += 1
        candidates.append({'path': d, 'ext': ext, 'kind': 'vendor_dir', 'inner_files': vendor_dirs[d]})
    for p in flat:
        if is_sidecar(p, cfg):
            census['%s (sidecar)' % extension(re.sub(r'\.[^.]+$', '', p), cfg)] += 1
            continue
        ext = extension(p, cfg)
        census[ext] += 1
        if ext in cfg.searchable_formats or ext in cfg.convertible_formats:
            candidates.append({'path': p, 'ext': ext, 'kind': 'file'})

    # One injection can be published several times: the vendor original plus the converted copy,
    # and sometimes under both peak/ and ccms_peak/. Those are representations of ONE acquisition,
    # so they collapse to a single row; the searchable copy is preferred as the representative.
    order = list(cfg.searchable_formats) + list(cfg.convertible_formats)

    def rank(a):
        top = a['path'].split('/')[0]
        pref = cfg.collection_preference
        return (order.index(a['ext']) if a['ext'] in order else len(order),
                pref.index(top) if top in pref else len(pref), len(a['path']))

    best = {}
    for a in sorted(candidates, key=rank):
        key = strip_wrappers(a['path'], cfg).rsplit('.', 1)[0]
        if key in best:
            best[key]['representations'] += 1
            best[key].setdefault('other_paths', []).append(a['path'])
            if a['ext'] not in best[key]['formats']:
                best[key]['formats'].append(a['ext'])
        else:
            a['representations'] = 1
            a['formats'] = [a['ext']]
            best[key] = a
    return list(best.values()), census


def scope_of(ext, cfg: Config):
    if ext in cfg.searchable_formats:
        return 'in_scope_searchable'
    if ext in cfg.convertible_formats:
        return 'in_scope_convert'
    return 'out_of_scope'


def converter_for(ext, vendor):
    return CONVERTERS.get((ext, (vendor or '').lower()), CONVERTERS.get((ext, ''), ''))


# -------------------------------------------------------------------- assay ----


def classify_assay(proxi, census, paths, charge_probe=None):
    """Metabolomics vs proteomics. Returns (label, confidence, evidence list).

    Ladder, in order of decisiveness; measured against known deposits of both kinds.
    """
    ev = []
    kws = [k.get('value', '') for k in (proxi.get('keywords') or [])]
    declared = [k.split(':', 1)[1].strip() for k in kws if k.startswith('DatasetType:')]
    if declared:
        ev.append('PROXI DatasetType keyword(s): %s' % declared)
        labs = {d.lower() for d in declared}
        if 'metabolomics' in labs and 'proteomics' in labs:
            return 'ambiguous', 'declared_both', ev
        if 'metabolomics' in labs:
            return 'metabolomics', 'declared', ev
        if 'proteomics' in labs:
            return 'proteomics', 'declared', ev
        return 'other', 'declared', ev

    lowered = [p.lower() for p in paths]
    hits = {e: sum(1 for p in lowered if p.endswith(e)) for e in PROTEOMICS_ONLY}
    hits = {k: v for k, v in hits.items() if v}
    if hits:
        ev.append('proteomics-only file types present: %s' % hits)
        return 'proteomics', 'search_files', ev
    ev.append('no DatasetType keyword and no proteomics-only file types')

    mods = [m.get('name', '') for m in (proxi.get('modifications') or [])]
    ptms = [m for m in mods if m and 'No PTMs' not in m]
    if ptms:
        ev.append('PROXI modifications list PTMs: %s' % ptms[:4])

    if charge_probe:
        resolves, ms2, det = charge_probe()
        charges = [d.get('precursor_charge') for d in det if (d.get('precursor_mz') or 0) > 0]
        ev.append('charge probe on MS2 scans: %s' % charges)
        if charges:
            multi = sum(1 for c in charges if (c or 0) >= 2)
            if multi >= max(1, len(charges) // 2):
                return 'proteomics', 'precursor_charge', ev
            return 'metabolomics', 'precursor_charge', ev

    if ptms:
        return 'proteomics', 'ptm_list_only', ev
    title = ((proxi.get('title') or '') + ' ' + (proxi.get('summary') or '')).lower()
    prot_words = sum(w in title for w in ('proteom', 'peptide', 'trypsin', 'tmt', 'itraq', 'phospho'))
    metab_words = sum(w in title for w in ('metabolom', 'lipidom', 'natural product', 'metabolite',
                                           'untargeted', 'gnps'))
    ev.append('title/summary word counts: proteomics=%d metabolomics=%d' % (prot_words, metab_words))
    if prot_words and not metab_words:
        return 'proteomics', 'title_only', ev
    if metab_words and not prot_words:
        return 'metabolomics', 'title_only', ev
    return 'ambiguous', 'no_signal', ev
