"""Scoring and validation. The model proposes; everything here decides whether it holds.

Each function takes a *proposal* (a regex, a CV term, a join key) and returns a verdict with
the statistic that produced it. No function here accepts a filled cell from a model directly.
"""
from __future__ import annotations

import csv, io, json, os, re
from collections import Counter, defaultdict

# ------------------------------------------------------------------ profile ---


def profile(stems, n_examples=30):
    """A statistical profile plus a stratified sample - never a raw dump of all filenames."""
    toks = [re.split(r'[_\-. ]', s) for s in stems]
    positions = []
    for i in range(max((len(t) for t in toks), default=0)):
        vals = [t[i] for t in toks if len(t) > i]
        if not vals:
            continue
        positions.append({
            'position': i, 'n': len(vals), 'cardinality': len(set(vals)),
            'alphabet': ''.join(sorted({'d' if c.isdigit() else 'a' if c.isalpha() else 'x'
                                        for v in vals for c in v})),
            'examples': sorted(set(vals))[:8],
        })
    by_len = defaultdict(list)
    for s, t in zip(stems, toks):
        by_len[len(t)].append(s)
    ex, per = [], max(1, n_examples // max(1, len(by_len)))
    for k in sorted(by_len):
        ex += by_len[k][:per]
    return {'n_files': len(stems), 'token_count_distribution': dict(Counter(len(t) for t in toks)),
            'positions': positions, 'examples': ex[:n_examples]}


def subtrees(paths):
    t = Counter(p.rsplit('/', 1)[0] if '/' in p else '' for p in paths)
    return [{'subtree': k, 'files': v} for k, v in t.most_common()]


# ------------------------------------------------------------------ grammar ---


def score_grammar(stems, patterns):
    """Run proposed regexes over every stem. Returns coverage and per-family assignment."""
    try:
        compiled = {name: re.compile(p) for name, p in patterns.items()}
    except re.error as e:
        return {'ok': False, 'error': 'bad regex: %s' % e, 'coverage': 0.0}
    assign, unmatched = {}, []
    for s in stems:
        for name, rx in compiled.items():
            m = rx.match(s)
            if m:
                assign[s] = {'family': name, **m.groupdict()}
                break
        else:
            unmatched.append(s)
    cov = len(assign) / len(stems) if stems else 0.0
    fields = defaultdict(Counter)
    for a in assign.values():
        for k, v in a.items():
            if k != 'family' and v is not None:
                fields[k][v] += 1
    return {'ok': True, 'coverage': cov, 'n_matched': len(assign), 'n_total': len(stems),
            'unmatched_examples': unmatched[:20], 'n_unmatched': len(unmatched),
            'family_counts': dict(Counter(a['family'] for a in assign.values())),
            'field_cardinality': {k: len(v) for k, v in fields.items()},
            'assignment': assign}


def coverage_by_subtree(assignment, path_of, stems):
    out = defaultdict(lambda: [0, 0])
    for s in stems:
        sub = path_of[s].rsplit('/', 1)[0]
        out[sub][1] += 1
        if s in assignment:
            out[sub][0] += 1
    return {k: {'matched': v[0], 'total': v[1], 'coverage': v[0] / v[1] if v[1] else 0}
            for k, v in out.items()}


def field_roles(scored):
    """join_key_candidate / factor / technical / noise, from cardinality alone."""
    n = scored['n_matched'] or 1
    roles = {}
    for fieldname, card in scored['field_cardinality'].items():
        ratio = card / n
        if ratio > 0.9:
            role = 'join_key_candidate'
        elif card == 1:
            role = 'noise'
        elif card <= 12:
            role = 'factor'
        else:
            role = 'technical'
        roles[fieldname] = {'cardinality': card, 'ratio_to_files': round(ratio, 3), 'role': role}
    return roles


# --------------------------------------------------------------------- join ---

NORMALIZERS = {
    'exact': lambda s: s,
    'lower': lambda s: s.lower(),
    'strip_zeros': lambda s: re.sub(r'(?<=[A-Za-z])0+(\d)', r'\1', s),
    'delim_swap': lambda s: re.sub(r'[-. ]', '_', s),
    'alnum_only': lambda s: re.sub(r'[^A-Za-z0-9]', '', s).lower(),
}


def test_join(file_keys, table_rows, column, cfg):
    """Every (normalization) pair tested, with match rate and injectivity against thresholds."""
    results = []
    files = set(file_keys)
    for norm_name, fn in NORMALIZERS.items():
        tvals = [str(r.get(column, '')).strip() for r in table_rows]
        tvals = [v for v in tvals if v]
        tnorm = Counter(fn(v) for v in tvals)
        fnorm = {fn(f) for f in files}
        matched = sum(1 for f in fnorm if f in tnorm)
        rate = matched / len(fnorm) if fnorm else 0.0
        rows_per_file = (sum(tnorm[f] for f in fnorm if f in tnorm) / matched) if matched else 0.0
        ok = rate >= cfg.join_match_rate and rows_per_file <= cfg.join_injectivity
        results.append({'normalization': norm_name, 'column': column, 'match_rate': round(rate, 4),
                        'rows_per_file': round(rows_per_file, 4), 'accept': ok,
                        'threshold': {'match_rate': cfg.join_match_rate,
                                      'rows_per_file': cfg.join_injectivity}})
    results.sort(key=lambda r: (-r['match_rate'], r['rows_per_file']))
    return results


def find_id_columns(table_rows, min_cardinality=0.5):
    """ID-like columns: high cardinality relative to row count."""
    if not table_rows:
        return []
    out = []
    for c in table_rows[0].keys():
        vals = [str(r.get(c, '')).strip() for r in table_rows]
        vals = [v for v in vals if v]
        if not vals:
            continue
        card = len(set(vals)) / len(vals)
        if card >= min_cardinality:
            out.append({'column': c, 'cardinality_ratio': round(card, 3), 'n': len(vals)})
    return sorted(out, key=lambda x: -x['cardinality_ratio'])


# ----------------------------------------------------------------------- cv ---


class ReduCV:
    """Allowed terms per ReDU column, fetched once and cached on disk."""

    TERMS = 'https://redu.gnps2.org/attribute/%s/attributeterms'
    ATTRS = 'https://redu.gnps2.org/attributes'

    def __init__(self, http, cache_dir):
        self.http = http
        self.dir = os.path.join(cache_dir, 'cv')
        os.makedirs(self.dir, exist_ok=True)
        self._mem = {}

    def attributes(self):
        txt = self.http.cached(self.ATTRS, 'cv/_attributes.json')
        return [a['attributename'] for a in json.loads(txt)]

    def terms(self, attr):
        if attr not in self._mem:
            txt = self.http.cached(self.TERMS % attr, 'cv/%s.json' % attr)
            try:
                self._mem[attr] = {t['attributeterm'] for t in json.loads(txt)}
            except Exception:
                self._mem[attr] = set()
        return self._mem[attr]

    def validate(self, attr, value):
        if value in ('', None):
            return True, 'empty'
        allowed = self.terms(attr)
        if not allowed:
            return True, 'no CV published for this column'
        return (value in allowed), ('in CV' if value in allowed else 'NOT in CV')

    def validate_table(self, rows, columns=None):
        """Returns a list of violations; an empty list is the only acceptable result."""
        bad = []
        cols = columns or (list(rows[0].keys()) if rows else [])
        for c in cols:
            allowed = self.terms(c)
            if not allowed:
                continue
            for v in {r.get(c, '') for r in rows}:
                if v and v not in allowed:
                    bad.append({'column': c, 'value': v})
        return bad

    # NOTE: there is deliberately no suggest()/shortlist helper here.
    # Choosing a controlled term is a reasoning task over evidence from several sources at once
    # (title, summary, keywords, species block, headers, filenames, submitter tables). Keyword
    # overlap looked helpful and was actively harmful: the phrase "resuspension in 50:50
    # MeOH:water" ranked water_ground and water_seawater top for a human stool deposit. The model
    # receives the whole vocabulary and reasons; code's only job is to check membership.


# ------------------------------------------------------------------ headers ---

HDR_PATTERNS = {
    'startTimeStamp': r'startTimeStamp="([^"]+)"',
    'instrument_model': r'accession="MS:100[0-9]{4}"\s+name="((?:Orbitrap|Exploris|Q Exactive|LTQ|timsTOF|Synapt|maXis|impact|TripleTOF)[^"]*)"',
    'instrument_serial': r'name="instrument serial number"[^>]*value="([^"]*)"',
    'ionization': r'name="(electrospray ionization|nanoelectrospray|atmospheric pressure chemical ionization|electron ionization)"',
    'analyzer': r'name="(orbitrap|quadrupole|time-of-flight|ion trap)"',
}


def parse_header(blob: bytes):
    t = blob.decode('utf-8', 'ignore')
    out = {k: (re.search(p, t).group(1) if re.search(p, t) else None) for k, p in HDR_PATTERNS.items()}
    pol = sorted(set(re.findall(r'name="(positive scan|negative scan)"', t)))
    out['polarities'] = pol
    out['polarity'] = pol[0] if len(pol) == 1 else ('mixed' if len(pol) > 1 else None)
    return out


def sample_per_subtree(paths, k):
    by = defaultdict(list)
    for p in paths:
        by[p.rsplit('/', 1)[0]].append(p)
    out = []
    for sub, ps in by.items():
        ps = sorted(ps)
        out += ps[:1] + (ps[len(ps) // 2:len(ps) // 2 + k - 1] if k > 1 else [])
    return out
