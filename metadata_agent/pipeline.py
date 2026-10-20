"""The orchestrator. Python calls every tool; the model is consulted at five points only.

Decision points (everything else is deterministic):
  D-ASSAY     only when the metabolomics/proteomics ladder is ambiguous
  D-GRAM      propose filename grammars; code scores coverage and may reject
  D-CV        map a quoted phrase to a controlled term; code validates membership
  D-CONFLICT  adjudicate a contradiction between two sources; code supplies both
  D-LINK      judge publication candidates; code enforces the hard-evidence rule
"""
from __future__ import annotations

import csv, io, json, os, re
from collections import Counter, defaultdict

from . import analysis as A
from . import core as C
from .engine import Decision, DecisionLog, build_engine

# Columns whose correct value can differ between the families inside one deposit: a solvent blank,
# a standards mix and a study sample are not the same kind of thing. Asking once per deposit is
# how blanks silently inherit study-sample metadata.
PER_FAMILY_COLUMNS = ['SampleType', 'SampleTypeSub1', 'NCBITaxonomy', 'UBERONBodyPartName',
                      'SampleCollectionMethod']
# Fallback only; the live schema is read from ReDU at run time.
FALLBACK_COLUMNS = [
    'ATTRIBUTE_DatasetAccession', 'filename', 'SampleType', 'SampleTypeSub1', 'NCBITaxonomy',
    'YearOfAnalysis', 'SampleCollectionMethod', 'SampleExtractionMethod', 'InternalStandardsUsed',
    'MassSpectrometer', 'IonizationSourceAndPolarity', 'ChromatographyAndPhase',
    'SubjectIdentifierAsRecorded', 'AgeInYears', 'BiologicalSex', 'UBERONBodyPartName',
    'TermsofPosition', 'HealthStatus', 'SmokingStatus', 'BodyMassIndex', 'Diet', 'DOIDCommonName',
    'ComorbidityListDOIDIndex', 'SampleCollectionDateandTime', 'Country', 'HumanPopulationDensity',
    'LatitudeandLongitude', 'DepthorAltitudeMeters', 'qiita_sample_name', 'UniqueSubjectID',
    'LifeStage', 'UBERONOntologyIndex', 'DOIDOntologyIndex', 'NCBIRank', 'NCBIDivision',
    'ENVOEnvironmentBiome', 'ENVOEnvironmentMaterial', 'ENVOBroadScale', 'ENVOLocalScale',
    'ENVOMediumScale', 'DataSource', 'USI',
]
# Two non-ReDU columns carried alongside the schema, appended after it so a strict ReDU
# consumer can still read the first N columns unchanged.
EXTRA_COLUMNS = ['file_format', 'indexed']
FORMAT_NAMES = {'.mzml': 'mzML', '.mzxml': 'mzXML', '.mgf': 'mgf', '.raw': 'raw', '.d': 'd',
                '.wiff': 'wiff', '.mzdata.xml': 'mzData'}

CV_COLUMNS = ['SampleType', 'SampleTypeSub1', 'SampleExtractionMethod', 'SampleCollectionMethod',
              'IonizationSourceAndPolarity', 'ChromatographyAndPhase', 'MassSpectrometer',
              'UBERONBodyPartName', 'LifeStage', 'HealthStatus', 'BiologicalSex']
TABLE_EXT = ('.csv', '.tsv', '.txt', '.xlsx')


class Run:
    def __init__(self, acc, outdir, cfg=None, engine='auto', summary_path=None,
                 extras=False, keep_log=False):
        self.acc, self.outdir = acc, outdir
        self.cfg = cfg or C.Config()
        os.makedirs(outdir, exist_ok=True)
        self.cfg.cache_dir = os.path.join(outdir, '_cache')
        self.http = C.Http(self.cfg)
        self.massive, self.cache = C.MassIVE(self.http), C.DatasetCache(self.http)
        self.usi, self.cv = C.UsiResolver(self.http), A.ReduCV(self.http, self.cfg.cache_dir)
        self.redu, self.ols = C.ReDU(self.http), C.OLS4(self.http)
        self.columns = list(FALLBACK_COLUMNS)
        self.engine = build_engine(engine, outdir) if isinstance(engine, str) else engine
        self.log = DecisionLog(os.path.join(outdir, 'decisions_%s.jsonl' % acc))
        self.summary_path, self.extras, self.keep_log = summary_path, extras, keep_log
        self.state = {'accession': acc, 'engine': self.engine.name}

    # ------------------------------------------------------------- stage 0 ---
    def acquire(self):
        try:
            self.columns = self.redu.schema()
            self.log.record('D-ACQ', 0, 'Which ReDU columns must this run account for?',
                            '%d columns read live from the ReDU dump header' % len(self.columns),
                            'The schema is read from ReDU with a Range request, not hardcoded, so a '
                            'column added upstream cannot silently go unfilled and unabstained.',
                            ['first/last: %s ... %s' % (self.columns[0], self.columns[-1])],
                            'A', ['schema'])
        except Exception as e:
            self.log.record('D-ACQ', 0, 'Which ReDU columns must this run account for?',
                            '%d columns from the built-in fallback list' % len(self.columns),
                            'Live schema fetch failed (%s); the fallback may be stale.' % str(e)[:120],
                            [], 'C', ['schema'])
        proxi = self.massive.proxi(self.acc)
        files = self.cache.filenames(self.acc)
        umri = {r['filepath']: r for r in self.cache.uniquemri(self.acc)}
        paths = [r['filepath'] for r in files]
        acqs, census = C.acquisitions(paths, self.cfg)
        in_scope = [a for a in acqs if C.scope_of(a['ext'], self.cfg) != 'out_of_scope']
        self.log.record('D-ACQ', 0, 'Which sources returned data?',
                        'PROXI + dataset cache (filename, uniquemri)',
                        'All three returned; %d cache rows, %d uniquemri rows.' % (len(files), len(umri)),
                        ['PROXI title: %s' % (proxi.get('title') or '')[:80],
                         'cache rows %d' % len(files), 'uniquemri rows %d' % len(umri)],
                        'A', ['all'])
        oos = {e: n for e, n in census.items() if C.scope_of(e, self.cfg) == 'out_of_scope'}
        self.log.record('D-FORMAT', 0, 'Format census after wrapper normalisation and vendor collapse',
                        '%d acquisitions in scope; %d out-of-scope files' % (len(in_scope), sum(oos.values())),
                        'Only %s are searchable and %s are convertible; everything else is counted, not rowed.'
                        % (list(self.cfg.searchable_formats), list(self.cfg.convertible_formats)),
                        ['census: %s' % dict(census), 'cache rows %d' % len(files)],
                        'A', ['row count', 'file_format'])
        for ext, n in sorted(oos.items()):
            self.log.abstain('files with extension %r' % (ext or '(none)'), 'out_of_scope_format',
                             '%d such files exist in the deposit and are deliberately not rowed.' % n,
                             stage=0)
        # spectrum counts
        counted = sum(1 for a in in_scope
                      if int(umri.get(a['path'], {}).get('spectra_ms2') or 0) > 0)
        self.state.update(proxi=proxi, files=files, umri=umri, paths=paths, acqs=acqs,
                          census=dict(census), in_scope=in_scope, ms2_counted=counted)
        self.log.record('D-ACQ', 0, 'Has GNPS2 parsed this deposit?',
                        '%d of %d in-scope acquisitions carry MS2 counts' % (counted, len(in_scope)),
                        'Zero counts mean not-yet-indexed, which is not the same as no MS/MS.',
                        ['uniquemri rows %d' % len(umri)], 'A', ['gnps2_ms2_spectra', 'searchable_now'])
        return self.state

    def classify_assay(self):
        s = self.state
        probe = None
        if s['in_scope']:
            first = s['in_scope'][0]['path']
            probe = lambda: self.usi.probe('mzspec:%s:%s' % (self.acc, first))
        label, conf, ev = C.classify_assay(s['proxi'], s['census'], s['paths'], probe)
        if label == 'ambiguous':
            d = Decision('D-ASSAY', 'Is this deposit metabolomics or proteomics?',
                         {'title': s['proxi'].get('title'), 'summary': s['proxi'].get('summary'),
                          'keywords': [k.get('value') for k in (s['proxi'].get('keywords') or [])],
                          'modifications': [m.get('name') for m in (s['proxi'].get('modifications') or [])],
                          'format_census': s['census'], 'ladder_evidence': ev},
                         options=['metabolomics', 'proteomics', 'both', 'unclear'])
            ans, tr = self.engine.ask(d)
            if ans:
                label, conf = ans['chosen'], 'model'
                self.log.record('D-ASSAY', 0, d.question, label, ans.get('basis', ''), ev, 'C',
                                ['whole run'], ans.get('options_considered'), ans.get('unresolved'),
                                decided_by=self.engine.name, transcript=tr)
            else:
                self.log.record('D-ASSAY', 0, d.question, 'ambiguous - not resolved',
                                'The deterministic ladder was inconclusive and no engine answer was available.',
                                ev, 'C', ['whole run'], decided_by='abstain', transcript=tr)
        else:
            self.log.record('D-ASSAY', 0, 'Is this deposit metabolomics or proteomics?', label,
                            'Resolved deterministically by the %s rung of the ladder.' % conf,
                            ev, 'A', ['whole run'])
        self.state['assay'] = label
        return label

    def sample_headers(self, k=None):
        s = self.state
        searchable = [a['path'] for a in s['in_scope']
                      if C.scope_of(a['ext'], self.cfg) == 'in_scope_searchable']
        picks = A.sample_per_subtree(searchable, k or self.cfg.header_sample_per_subtree)
        heads = {}
        for p in picks:
            try:
                heads[p] = A.parse_header(self.massive.fetch_bytes(self.acc, p))
            except Exception as e:
                heads[p] = {'error': str(e)[:120]}
        models = Counter(h.get('instrument_model') for h in heads.values() if h.get('instrument_model'))
        pols = Counter(h.get('polarity') for h in heads.values() if h.get('polarity'))
        ions = Counter(h.get('ionization') for h in heads.values() if h.get('ionization'))
        self.state.update(headers=heads, header_models=models, header_pols=pols, header_ions=ions)
        self.log.record('D-ACQ', 0, 'What do the raw-file headers say?',
                        'model=%s polarity=%s ionization=%s' % (dict(models), dict(pols), dict(ions)),
                        'Sampled %d files across %d subtrees; headers outrank prose for these columns.'
                        % (len(heads), len({p.rsplit('/', 1)[0] for p in picks})),
                        ['sampled: %s' % picks[:6]], 'A',
                        ['MassSpectrometer', 'IonizationSourceAndPolarity', 'YearOfAnalysis'])
        if len(pols) > 1:
            self.log.record('D-CONFLICT', 3, 'Polarity differs across sampled files',
                            'mixed polarity recorded', 'Headers disagree; this is a finding, not noise.',
                            ['polarities: %s' % dict(pols)], 'A', ['IonizationSourceAndPolarity'])
        return heads

    # ------------------------------------------------------------- stage 1 ---
    def induce_grammar(self):
        s = self.state
        searchable = [a for a in s['in_scope']]
        stems = [C.strip_wrappers(a['path'], self.cfg).rsplit('.', 1)[0] for a in searchable]
        path_of = {C.strip_wrappers(a['path'], self.cfg).rsplit('.', 1)[0]: a['path'] for a in searchable}
        prof = A.profile(stems)
        prof['subtrees'] = A.subtrees([a['path'] for a in searchable])[:10]
        best, tr_all = None, []
        d = Decision('D-GRAM',
                     'Propose named regex families that together explain every acquisition filename.',
                     {'profile': prof,
                      'instruction': 'Return {"chosen": {"<family>": "<python regex with named groups>"}, '
                                     '"basis": "..."}. Anchor each regex with ^ and $. Use named groups '
                                     'for sample ids, plate, well, replicate, injection index, QC role. '
                                     'Several families are normal: study samples, QC, blanks.'},
                     required_keys=('chosen', 'basis'))

        def validator(ans):
            pats = ans.get('chosen')
            if not isinstance(pats, dict) or not pats:
                return False, 'chosen must be a non-empty object of {family: regex}'
            scored = A.score_grammar(stems, pats)
            if not scored['ok']:
                return False, scored['error']
            ans['_scored'] = scored
            if scored['coverage'] < self.cfg.grammar_coverage:
                return False, ('coverage %.3f is below the %.2f threshold; unmatched examples: %s'
                               % (scored['coverage'], self.cfg.grammar_coverage,
                                  scored['unmatched_examples'][:10]))
            return True, 'coverage %.3f' % scored['coverage']

        ans, tr = self.engine.ask(d, validator)
        tr_all += tr
        if ans:
            scored = ans['_scored']
            roles = A.field_roles(scored)
            self.log.record('D-GRAM', 1, d.question, json.dumps(ans['chosen']), ans.get('basis', ''),
                            ['coverage %.3f over %d files' % (scored['coverage'], scored['n_total']),
                             'families: %s' % scored['family_counts'],
                             'per-subtree coverage: %s' % A.coverage_by_subtree(scored['assignment'], path_of, stems)],
                            'C', list(roles), ans.get('options_considered'), ans.get('unresolved'),
                            decided_by=self.engine.name, transcript=tr_all)
            self.log.record('D-FIELD', 1, 'Role of each induced field', json.dumps(roles),
                            'Assigned from cardinality relative to file count.',
                            ['%s' % roles], 'C', list(roles))
            best = {'patterns': ans['chosen'], 'scored': scored, 'roles': roles}
        else:
            self.log.record('D-GRAM', 1, d.question, 'no grammar accepted',
                            'No proposal reached the %.0f%% coverage threshold, or no engine was available. '
                            'Structural variables stay empty rather than being guessed.'
                            % (100 * self.cfg.grammar_coverage), [], 'A',
                            ['plate_id', 'well_position', 'injection_order', 'qc_role'],
                            decided_by='abstain', transcript=tr_all)
            for col in ('plate_id', 'well_position', 'replicate_index', 'qc_role'):
                self.log.abstain(col, 'no_source', 'No filename grammar cleared the coverage threshold.',
                                 stage=1)
        self.state['grammar'] = best
        self.state['stems'] = stems
        self.state['path_of'] = path_of
        return best

    # ------------------------------------------------------------- stage 1b --
    def try_external_tables(self):
        """Submitter-supplied tables inside the deposit, joined under the hard thresholds."""
        s = self.state
        tables = [p for p in s['paths'] if p.lower().endswith(TABLE_EXT)
                  and not p.lower().endswith(('.mzml.scans', '.mzxml.scans'))]
        joined = []
        for t in tables[:4]:
            if t.lower().endswith('.xlsx'):
                continue
            try:
                blob = self.massive.fetch_bytes(self.acc, t, nbytes=2_000_000)
            except Exception as e:
                self.log.record('D-JOIN', 1, 'Fetch %s' % t, 'failed', str(e)[:150], [], 'A', [])
                continue
            text = blob.decode('utf-8', 'ignore')
            rows, header_row = self._sniff(text)
            if not rows:
                continue
            idcols = A.find_id_columns(rows)
            best = None
            for c in idcols[:4]:
                res = A.test_join(s['stems'], rows, c['column'], self.cfg)
                top = res[0]
                self.log.record('D-JOIN', 1, 'Join %s on column %r' % (t, c['column']),
                                'accept' if top['accept'] else 'reject',
                                'Best normalization %r gave match rate %.3f with %.3f rows per file '
                                '(thresholds %.2f / %.2f).'
                                % (top['normalization'], top['match_rate'], top['rows_per_file'],
                                   self.cfg.join_match_rate, self.cfg.join_injectivity),
                                ['all normalizations: %s' % json.dumps(res)], 'A', ['external metadata'])
                if top['accept'] and (best is None or top['match_rate'] > best[0]['match_rate']):
                    best = (top, c['column'], rows, header_row)
            if best:
                joined.append({'table': t, 'column': best[1], 'result': best[0],
                               'rows': best[2], 'header_offset': best[3]})
        self.state['external'] = joined
        return joined

    @staticmethod
    def _sniff(text):
        """Instrument sequence files and plate maps often carry a preamble line."""
        lines = [l for l in text.splitlines() if l.strip()]
        for skip in (0, 1, 2):
            try:
                rows = list(csv.DictReader(io.StringIO('\n'.join(lines[skip:]))))
            except Exception:
                continue
            if rows and len(rows[0]) > 2 and sum(1 for k in rows[0] if k and k.strip()) > 2:
                return rows, skip
        return [], 0

    # ------------------------------------------------------------- stage 2 ---
    def evidence_bundle(self):
        """Everything known about the deposit, from every source, in one object.

        Controlled-term selection reasons over all of it at once - the instrument in a header can
        settle a chromatography question that the prose never mentions, and a submitter table can
        contradict the title. Handing the model one source at a time would hide exactly that.
        """
        s = self.state
        heads = s.get('headers') or {}
        ext = []
        for e in s.get('external', []):
            rows = e.get('rows') or []
            ext.append({'table': e['table'], 'joined_on': e['column'],
                        'match_rate': e['result']['match_rate'],
                        'columns': list(rows[0].keys()) if rows else [],
                        'sample_rows': rows[:3]})
        return {
            'repository_record': {
                'title': s['proxi'].get('title'), 'summary': s['proxi'].get('summary'),
                'keywords': [k.get('value') for k in (s['proxi'].get('keywords') or [])],
                'species': s['proxi'].get('species'),
                'instruments': [i.get('name') for i in (s['proxi'].get('instruments') or [])],
                'publications': [p.get('name') for p in (s['proxi'].get('publications') or [])]},
            'raw_file_headers': {
                'instrument_models': dict(s.get('header_models') or {}),
                'polarities': dict(s.get('header_pols') or {}),
                'ionization_sources': dict(s.get('header_ions') or {}),
                'analyzers': dict(Counter(h.get('analyzer') for h in heads.values() if h.get('analyzer'))),
                'acquisition_timestamps': sorted({h.get('startTimeStamp') for h in heads.values()
                                                  if h.get('startTimeStamp')})[:6],
                'n_sampled': len(heads)},
            'file_layout': {'format_census': s.get('census'),
                            'subtrees': A.subtrees([a['path'] for a in s['in_scope']])[:8],
                            'example_filenames': [a['path'].rsplit('/', 1)[-1] for a in s['in_scope'][:12]]},
            'submitter_tables': ext,
            'grammar': (s.get('grammar') or {}).get('patterns'),
            'note': 'Headers outrank prose for instrument, polarity and acquisition date. '
                    'A value stated nowhere in these sources must be abstained, not inferred.',
        }

    def ground_to_cv(self):
        """Every controlled-term choice is a model decision over the full vocabulary.

        Code does not pre-filter, rank or shortlist the vocabulary: that is reasoning, and keyword
        matching gets it wrong in exactly the cases that matter. Code supplies the evidence, lists
        the permitted terms, and then checks that whatever comes back is one of them.
        """
        s = self.state
        bundle = self.evidence_bundle()
        filled = {}
        for col in ['MassSpectrometer', 'IonizationSourceAndPolarity', 'ChromatographyAndPhase',
                    'SampleType', 'SampleTypeSub1', 'UBERONBodyPartName', 'SampleExtractionMethod',
                    'SampleCollectionMethod']:
            allowed = sorted(self.cv.terms(col))
            if not allowed:
                continue
            d = Decision('D-CV',
                         'Which %s term do these sources support? Abstain unless the evidence states it.' % col,
                         {'column': col, 'evidence': bundle,
                          'allowed_terms': allowed, 'n_allowed_terms': len(allowed),
                          'instruction': 'Return {"chosen": "<one exact term from allowed_terms, or ABSTAIN>", '
                                         '"basis": "...", "quote": "<verbatim supporting phrase or header value>", '
                                         '"source": "<which evidence block decided it>"}. '
                                         'Weigh all evidence blocks together. ABSTAIN is the correct answer '
                                         'whenever the sources do not settle it.'},
                         required_keys=('chosen', 'basis'))

            def cv_validator(ans, _col=col):
                if ans.get('chosen') == 'ABSTAIN':
                    return True, 'abstained'
                ok, why = self.cv.validate(_col, ans.get('chosen'))
                return ok, '%r %s' % (ans.get('chosen'), why)

            ans, tr = self.engine.ask(d, cv_validator)
            if not ans:
                self.log.abstain(col, 'no_source',
                                 'No engine answer was available, so the column is left empty rather '
                                 'than filled by a coded guess.', stage=2)
                continue
            if ans['chosen'] == 'ABSTAIN':
                self.log.record('D-ABSTAIN', 2, 'Why is %s empty?' % col, 'empty (no_source)',
                                ans.get('basis', 'The evidence does not state it.'),
                                ['quote: %s' % ans.get('quote', '')], 'A', [col],
                                decided_by=self.engine.name, transcript=tr)
                continue
            filled[col] = ans['chosen']
            self.log.record('D-CV', 2, d.question, ans['chosen'], ans.get('basis', ''),
                            ['quote: %s' % ans.get('quote', ''), 'source: %s' % ans.get('source', '')],
                            'B', [col], ans.get('options_considered'), ans.get('unresolved'),
                            decided_by=self.engine.name, transcript=tr)

        # Structured, non-vocabulary facts stay in code: these are reads, not judgements.
        years = {h.get('startTimeStamp', '')[:4] for h in (s.get('headers') or {}).values()
                 if h.get('startTimeStamp')}
        if len(years) == 1:
            filled['YearOfAnalysis'] = years.pop()
            self.log.record('D-CV', 2, 'Fill YearOfAnalysis', filled['YearOfAnalysis'],
                            'Read from mzML startTimeStamp; a date, not a vocabulary term.',
                            ['sampled headers agree on one year'], 'A', ['YearOfAnalysis'])
        # Taxonomy comes from the PROXI species block, a structured field, never from prose.
        # Per-family overrides (a blank contains no organism) are applied later.
        tax = {}
        for block in (s['proxi'].get('species') or []):
            for item in block:
                if item.get('accession') == 'MS:1001467':
                    tax['taxid'] = item.get('value')
                if item.get('accession') == 'MS:1001469':
                    tax['name'] = item.get('value')
        if tax.get('taxid') and tax.get('name') and 'NCBITaxonomy' in self.columns:
            filled['NCBITaxonomy'] = '%s|%s' % (tax['taxid'], tax['name'])
            self.log.record('D-CV', 2, 'Fill NCBITaxonomy', filled['NCBITaxonomy'],
                            'Read from the PROXI species block.', ['PROXI species: %s' % tax], 'A',
                            ['NCBITaxonomy'])
        self.state['dataset_level'] = filled
        return filled

    # ------------------------------------------------------- stage 2, per family
    def _external_lookup(self):
        """stem -> row of whichever submitter table passed the join thresholds."""
        out = {}
        for e in self.state.get('external', []):
            fn = A.NORMALIZERS[e['result']['normalization']]
            idx = {}
            for r in e['rows']:
                v = str(r.get(e['column'], '')).strip()
                if v:
                    idx[fn(v)] = r
            for stem in self.state['stems']:
                if fn(stem) in idx:
                    out[stem] = idx[fn(stem)]
        return out

    def per_family_cv(self):
        """Ask once per grammar family, not once per deposit.

        A solvent blank, a standards mix and a study sample are different materials; a single
        dataset-level answer is how blanks inherit study-sample metadata.
        """
        s = self.state
        g = s.get('grammar')
        if not g:
            self.log.record('D-CV', 2, 'Per-family controlled terms', 'skipped',
                            'No grammar was accepted, so there are no families to differentiate.',
                            [], 'A', PER_FAMILY_COLUMNS)
            self.state['family_values'] = {}
            return {}
        assignment = g['scored']['assignment']
        fams = defaultdict(list)
        for stem, a in assignment.items():
            fams[a['family']].append(stem)
        ext = self._external_lookup()
        size_of = {r['filepath']: int(r.get('size') or 0) for r in s['files']}
        cols = [c for c in PER_FAMILY_COLUMNS if c in self.columns]
        vocab = {c: sorted(self.cv.terms(c)) for c in cols}
        vocab = {c: v for c, v in vocab.items() if v}
        bundle = self.evidence_bundle()

        # A grammar family can still mix materials: blanks named alike may be typed differently in
        # the submitter's own sequence table. Where a joined table splits a family into a handful of
        # roles, ask per role instead of per family - otherwise a QC-typed blank inherits the
        # answer given for ordinary blanks.
        groups = {}
        for fam, stems in sorted(fams.items()):
            split_col, best = None, None
            cols_in_table = set()
            for st in stems:
                cols_in_table |= set(ext.get(st, {}).keys())
            for col in sorted(c for c in cols_in_table if c and c.strip()):
                vals = {str(ext[st].get(col, '')).strip() for st in stems if st in ext}
                vals.discard('')
                if 2 <= len(vals) <= 5 and (best is None or len(vals) < best):
                    split_col, best = col, len(vals)
            if split_col:
                for st in stems:
                    v = str(ext.get(st, {}).get(split_col, '')).strip() or '(absent)'
                    groups.setdefault((fam, split_col, v), []).append(st)
            else:
                groups.setdefault((fam, None, None), []).append(stems[0]) if False else None
                groups[(fam, None, None)] = stems
        if any(k[1] for k in groups):
            self.log.record('D-FIELD', 2, 'Does a submitter table subdivide any grammar family?',
                            '; '.join('%s/%s=%s (%d files)' % (f, c, v, len(st))
                                      for (f, c, v), st in sorted(groups.items()) if c),
                            'A joined table splits at least one family into roles, so controlled terms '
                            'are chosen per role rather than per family.',
                            ['groups: %s' % {('%s|%s=%s' % k if k[1] else k[0]): len(v)
                                             for k, v in sorted(groups.items())}],
                            'A', PER_FAMILY_COLUMNS)

        out, stem_group = {}, {}
        for (fam, split_col, split_val), stems in sorted(groups.items()):
            group_key = fam if not split_col else '%s|%s=%s' % (fam, split_col, split_val)
            stem_group.update({st: group_key for st in stems})
            sample = sorted(stems)[:10]
            sizes = sorted(round(size_of.get(s['path_of'][st], 0) / 1e6, 1)
                           for st in stems if st in s['path_of'])
            d = Decision('D-CV',
                         'For the %r group (%d files), which values do the sources support? '
                         'Abstain per column unless the evidence states it.' % (group_key, len(stems)),
                         {'group': group_key, 'grammar_family': fam,
                          'split_by': ({'column': split_col, 'value': split_val} if split_col else None),
                          'n_files': len(stems), 'example_filenames': sample,
                          'file_size_MB': {'min': sizes[0] if sizes else None,
                                           'median': sizes[len(sizes) // 2] if sizes else None,
                                           'max': sizes[-1] if sizes else None},
                          'all_groups': {('%s|%s=%s' % k if k[1] else k[0]): len(v) for k, v in sorted(groups.items())},
                          'submitter_table_rows_for_this_family':
                              [ext[st] for st in sample if st in ext][:4],
                          'dataset_level_values': s.get('dataset_level', {}),
                          'columns_to_fill': cols, 'allowed_terms': vocab,
                          'evidence': bundle,
                          'instruction': 'Return {"chosen": {"<column>": "<exact allowed term, or '
                                         'ABSTAIN>"}, "basis": "...", "quote": "..."}. One entry per '
                                         'column in columns_to_fill. NCBITaxonomy is free text of the '
                                         'form "taxid|Scientific name" - use ABSTAIN for material that '
                                         'contains no organism, such as a solvent blank or a chemical '
                                         'standards mix.'},
                         required_keys=('chosen', 'basis'))

            def validator(ans, _cols=cols, _vocab=vocab):
                ch = ans.get('chosen')
                if not isinstance(ch, dict):
                    return False, 'chosen must be an object of {column: term}'
                bad = [k for k in ch if k not in _cols]
                if bad:
                    return False, 'unknown columns %s' % bad
                for k, v in ch.items():
                    if v == 'ABSTAIN' or v == '':
                        continue
                    if k in _vocab and v not in _vocab[k]:
                        return False, '%r is not an allowed %s term' % (v, k)
                    if k == 'NCBITaxonomy' and not re.match(r'^\d+\|\S', str(v)):
                        return False, 'NCBITaxonomy must look like "9606|Homo sapiens"'
                return True, 'ok'

            ans, tr = self.engine.ask(d, validator)
            if not ans:
                self.log.record('D-ABSTAIN', 2, 'Per-group values for %r' % group_key, 'empty (no_source)',
                                'No engine answer; this family keeps only dataset-level values.',
                                ['%d files' % len(stems)], 'A', cols,
                                decided_by='abstain', transcript=tr)
                continue
            vals = {k: v for k, v in ans['chosen'].items() if v and v != 'ABSTAIN'}
            out[group_key] = vals
            self.log.record('D-CV', 2, d.question, json.dumps(vals), ans.get('basis', ''),
                            ['quote: %s' % ans.get('quote', ''), 'files: %d' % len(stems),
                             'size MB median: %s' % (sizes[len(sizes) // 2] if sizes else None)],
                            'B', cols, ans.get('options_considered'), ans.get('unresolved'),
                            decided_by=self.engine.name, transcript=tr)
            for c in cols:
                if c not in vals:
                    self.log.record('D-ABSTAIN', 2, 'Why is %s empty for %r?' % (c, group_key),
                                    'empty (no_source)', ans.get('basis', 'Not stated for this family.'),
                                    [], 'A', [c], decided_by=self.engine.name)
        self.state['family_values'] = out
        self.state['stem_group'] = stem_group
        return out

    # ------------------------------------------------- stage 2, ontology lookup
    def resolve_ontologies(self):
        """Ontology ids and taxonomic lineage come from OLS4, not from recall.

        The model is consulted only when OLS4 returns several candidates for a label.
        """
        s = self.state
        ds = s.get('dataset_level', {})
        fams = s.get('family_values', {})
        onto = {'by_uberon': {}, 'by_taxonomy': {}}

        labels = {v.get('UBERONBodyPartName') for v in list(fams.values()) + [ds]}
        for label in sorted(x for x in labels if x):
            chosen, cands = self.ols.resolve(label, 'uberon')
            if not chosen and cands:
                d = Decision('D-ONTO', 'Which OLS4 UBERON term is %r?' % label,
                             {'label': label, 'candidates': cands,
                              'context': {'title': s['proxi'].get('title'),
                                          'summary': s['proxi'].get('summary')},
                              'instruction': 'Return {"chosen": "<obo_id>", "basis": "..."} or '
                                             '{"chosen": "ABSTAIN", ...}'},
                             options=[c['obo_id'] for c in cands] + ['ABSTAIN'])
                ans, tr = self.engine.ask(d)
                if ans and ans['chosen'] != 'ABSTAIN':
                    chosen = next(c for c in cands if c['obo_id'] == ans['chosen'])
                    self.log.record('D-ONTO', 2, d.question, chosen['obo_id'], ans.get('basis', ''),
                                    ['OLS4 candidates: %s' % [c['obo_id'] for c in cands]], 'B',
                                    ['UBERONOntologyIndex'], decided_by=self.engine.name, transcript=tr)
            if chosen:
                onto['by_uberon'][label] = chosen['obo_id']
                self.log.record('D-ONTO', 2, 'OLS4 lookup for UBERON label %r' % label,
                                chosen['obo_id'],
                                'Resolved against the EBI Ontology Lookup Service, not from memory.',
                                ['OLS4 label: %s' % chosen['label'],
                                 'candidates returned: %d' % len(cands)], 'A', ['UBERONOntologyIndex'])
            else:
                self.log.abstain('UBERONOntologyIndex', 'ambiguous',
                                 'OLS4 returned %d candidates for %r and none was selected.'
                                 % (len(cands), label), stage=2)

        taxa = {v.get('NCBITaxonomy') for v in list(fams.values()) + [ds]}
        for tax in sorted(x for x in taxa if x):
            taxid = tax.split('|')[0].strip()
            if not taxid.isdigit():
                continue
            try:
                lineage, own_rank = self.ols.lineage(taxid)
            except Exception as e:
                self.log.abstain('NCBI lineage columns', 'no_source',
                                 'OLS4 lineage lookup failed for taxid %s (%s).' % (taxid, str(e)[:80]),
                                 stage=2)
                continue
            cols = {}
            for redu_col, rank in [('NCBISuperkingdom', 'superkingdom'), ('NCBIKingdom', 'kingdom'),
                                   ('NCBIPhylum', 'phylum'), ('NCBIClass', 'class'),
                                   ('NCBIOrder', 'order'), ('NCBIFamily', 'family'),
                                   ('NCBIGenus', 'genus'), ('NCBISpecies', 'species')]:
                if redu_col in self.columns and lineage.get(rank):
                    cols[redu_col] = lineage[rank]
            if own_rank and 'NCBIRank' in self.columns:
                cols['NCBIRank'] = own_rank
            onto['by_taxonomy'][tax] = cols
            self.log.record('D-ONTO', 2, 'OLS4 NCBITaxon lineage for %s' % tax, json.dumps(cols),
                            'Ancestors and ranks read from OLS4 (ncbitaxon + taxrank); nothing recalled.',
                            ['ranks returned: %s' % sorted(lineage)], 'A', sorted(cols))
            missing = [c for c in ('NCBIDivision',) if c in self.columns]
            for c in missing:
                self.log.abstain(c, 'no_source',
                                 'NCBI division is not a rank and OLS4 does not expose it; it would '
                                 'have to come from the NCBI taxonomy dump.', stage=2)
        self.state['ontology'] = onto
        return onto

    # ------------------------------------------------------------- stage 4 ---
    def gate(self):
        s = self.state
        empty_subject = [c for c in ('AgeInYears', 'BiologicalSex', 'HealthStatus', 'Diet',
                                     'SubjectIdentifierAsRecorded')
                         if c not in s.get('dataset_level', {})]
        roles = (s.get('grammar') or {}).get('roles') or {}
        has_key = any(r['role'] == 'join_key_candidate' for r in roles.values())
        pubs = [p.get('name', '') for p in (s['proxi'].get('publications') or [])]
        handle = bool((s['proxi'].get('title') or '').strip()) or any('no associated' not in p.lower() for p in pubs)
        conds = {'missing_are_subject_level': bool(empty_subject),
                 'join_key_candidate_exists': has_key,
                 'searchable_handle': handle}
        opened = all(conds.values())
        self.log.record('D-GATE', 4, 'Does this deposit earn a publication search?',
                        'OPEN' if opened else 'DECLINED',
                        'All three conditions must hold; this is computed, not judged.',
                        ['%s = %s' % (k, v) for k, v in conds.items()] + ['publications: %s' % pubs],
                        'A', ['stage 5 entry'])
        self.state['gate'] = {'opened': opened, 'conditions': conds, 'missing': empty_subject}
        return self.state['gate']

    def link_publication(self, candidates=None):
        """Stage 5. Candidates must be supplied by the caller (the package does no web search)."""
        if not self.state.get('gate', {}).get('opened'):
            return None
        if not candidates:
            self.log.record('D-LINK', 5, 'Is there a publication describing these samples?',
                            'no_confident_link',
                            'The gate opened but no candidate documents were supplied to this run, so '
                            'nothing could clear the hard-evidence bar.', [], 'A',
                            self.state['gate']['missing'])
            return None
        d = Decision('D-LINK', 'Which candidate, if any, clears the hard-evidence bar?',
                     {'accession': self.acc, 'title': self.state['proxi'].get('title'),
                      'candidates': candidates,
                      'rule': 'Accept only if the accession or dataset title appears in the paper, or '
                              'methods match on >=3 independent axes AND an author overlaps the '
                              'depositor list. Author overlap alone is never sufficient.'},
                     options=[c.get('id') for c in candidates] + ['no_confident_link'])
        ans, tr = self.engine.ask(d)
        chosen = ans['chosen'] if ans else 'no_confident_link'
        self.log.record('D-LINK', 5, d.question, chosen,
                        (ans or {}).get('basis', 'No engine answer; defaulted to no link.'),
                        ['candidates: %s' % [c.get('id') for c in candidates]], 'A',
                        self.state['gate']['missing'], (ans or {}).get('options_considered'),
                        decided_by=self.engine.name if ans else 'abstain', transcript=tr)
        return None if chosen == 'no_confident_link' else chosen

    # ------------------------------------------------------------- stage 7 ---
    def emit(self):
        s = self.state
        ds = s.get('dataset_level', {})
        gram = (s.get('grammar') or {}).get('scored', {}).get('assignment', {})
        rows, side = [], []
        for a in s['in_scope']:
            stem = C.strip_wrappers(a['path'], self.cfg).rsplit('.', 1)[0]
            u = s['umri'].get(a['path'], {})
            searchable = C.scope_of(a['ext'], self.cfg) == 'in_scope_searchable'
            g = gram.get(stem, {})
            rec = {c: '' for c in self.columns}
            rec.update({'ATTRIBUTE_DatasetAccession': self.acc, 'filename': 'f.' + a['path'],
                        'USI': 'mzspec:%s:%s' % (self.acc, a['path']), 'DataSource': 'GNPS'})
            rec.update({k: v for k, v in ds.items() if k in rec})
            # per-family values override dataset-level ones, and may also clear them
            fam_values = s.get('family_values') or {}
            group_key = (s.get('stem_group') or {}).get(stem)
            if fam_values and group_key in fam_values:
                for k in PER_FAMILY_COLUMNS:
                    if k in rec:
                        rec[k] = fam_values[group_key].get(k, '')
            # ontology expansions follow whichever identifier the row ended up with
            onto = s.get('ontology') or {}
            idx = onto.get('by_uberon', {}).get(rec.get('UBERONBodyPartName'))
            if idx and 'UBERONOntologyIndex' in rec:
                rec['UBERONOntologyIndex'] = idx
            for k, v in (onto.get('by_taxonomy', {}).get(rec.get('NCBITaxonomy')) or {}).items():
                if k in rec:
                    rec[k] = v
            ms1 = int(u.get('spectra_ms1') or 0)
            ms2 = int(u.get('spectra_ms2') or 0)
            fmts = a.get('formats') or [a['ext']]
            rec['file_format'] = '+'.join(FORMAT_NAMES.get(e, e.lstrip('.')) for e in fmts)
            # 'indexed' is about GNPS2 having parsed the file, not about the data being absent:
            # an unparsed deposit still holds MS/MS, it just has no counts yet.
            rec['indexed'] = 'yes' if ms2 > 0 else ('ms1_only' if ms1 > 0 else 'no')
            rows.append(rec)
            vendor = (u.get('instrument_vendor') or '')
            side.append({'filename': rec['filename'], 'stem': stem, 'file_format': a['ext'].lstrip('.'),
                         'searchable_now': str(searchable and int(u.get('spectra_ms2') or 0) > 0).lower(),
                         'conversion_required': str(not searchable).lower(),
                         'conversion_recipe': '' if searchable else C.converter_for(a['ext'], vendor),
                         'gnps2_ms1_spectra': u.get('spectra_ms1', ''),
                         'gnps2_ms2_spectra': u.get('spectra_ms2', ''),
                         'gnps2_classification': u.get('classification', ''),
                         'grammar_family': g.get('family', ''),
                         **{k: v for k, v in g.items() if k != 'family'},
                         'acquisition_datetime': (s.get('headers', {}).get(a['path'], {}) or {}).get('startTimeStamp', ''),
                         'polarity': (s.get('headers', {}).get(a['path'], {}) or {}).get('polarity', ''),
                         'inner_files': a.get('inner_files', '')})
        violations = self.cv.validate_table(rows, CV_COLUMNS)
        filled = [c for c in self.columns if any(r.get(c) for r in rows)]
        for c in self.columns:
            if c not in filled and c not in ('ATTRIBUTE_DatasetAccession',):
                self.log.abstain(c, 'no_source', 'Nothing in the repository record speaks to this column.')

        redu_path = os.path.join(self.outdir, 'redu_%s.csv' % self.acc)
        self._write_csv(redu_path, rows, self.columns + EXTRA_COLUMNS)
        summary = self.summary_row(rows, side, filled, violations)
        self._write_csv(os.path.join(self.outdir, 'decision_summary_%s.csv' % self.acc),
                        [summary], list(summary))
        if self.summary_path:
            self._upsert_csv(self.summary_path, summary, 'accession')
        if self.extras:
            self._write_csv(os.path.join(self.outdir, 'derived_%s.csv' % self.acc), side,
                            sorted({k for r in side for k in r}))
        return {'rows': len(rows), 'cv_violations': len(violations), 'decisions': self.log.n,
                'filled_columns': len(filled), 'outputs': [os.path.basename(redu_path),
                                                           'decision_summary_%s.csv' % self.acc]}

    # ------------------------------------------------------- the second doc ---
    def summary_row(self, rows, side, filled, violations):
        """One row per dataset: what was decided, on what information, in a few sentences.

        Composed from the decision log rather than written freely, so the prose cannot drift from
        what the run actually did.
        """
        s = self.state
        g = s.get('grammar') or {}
        cov = (g.get('scored') or {}).get('coverage')
        fams = (g.get('scored') or {}).get('family_counts') or {}
        groups = s.get('family_values') or {}
        ext = s.get('external') or []
        census = s.get('census') or {}
        in_scope_exts = sorted({a['ext'] for a in s.get('in_scope', [])})
        oos = {k: v for k, v in census.items() if k not in in_scope_exts}
        heads = s.get('headers') or {}
        gate = s.get('gate') or {}
        n_model = sum(1 for r in self._log_records() if r.get('decided_by') not in ('code', 'abstain'))
        n_abst = sum(1 for r in self._log_records() if r['id'].startswith('D-ABSTAIN'))
        conv = sum(1 for x in side if x.get('conversion_required') == 'true')
        searchable = sum(1 for x in side if x.get('searchable_now') == 'true')

        sources, missing = self._sources()

        assay_ev = next((r['evidence'][0] for r in self._log_records()
                         if r['id'].startswith('D-ASSAY') and r.get('evidence')), '')
        sent = []
        if not rows:
            sent.append('Classified as %s on %s, so metabolomics recovery was not run and no rows '
                        'were emitted.' % (s.get('assay'), assay_ev or 'the available evidence'))
            sent.append('The decision used only the %s; nothing was downloaded and the deposit %s.'
                        % (' and '.join(sources) if sources else 'repository record',
                           'is not held by the GNPS2 file cache' if not s.get('files')
                           else 'holds %d cache paths' % len(s['files'])))
            return dict(self._summary_fields(rows, side, filled, violations), how_decided=' '.join(sent))
        sent.append('Classified as %s on %s.' % (s.get('assay'), assay_ev or 'the available evidence'))
        part = ('%d acquisitions were rowed from %d cache paths after excluding %s'
                % (len(rows), len(s.get('files', [])),
                   ', '.join('%d %s' % (v, k) for k, v in sorted(oos.items())) or 'nothing'))
        if cov is not None:
            part += ('; a filename grammar matching %.0f%% of them split the deposit into %s'
                     % (100 * cov, ', '.join('%s (%d)' % (k, v) for k, v in sorted(fams.items()))))
        if ext:
            e = ext[0]
            part += ('; the submitter table %s joined on %r at %.1f%% with %.2f rows per file'
                     % (e['table'].rsplit('/', 1)[-1], e['column'],
                        100 * e['result']['match_rate'], e['result']['rows_per_file']))
        sent.append(part + '.')
        sent.append('%d of %d ReDU columns were filled and %d abstentions were recorded, with %d '
                    'controlled-term choice(s) made by the %s engine over the full vocabulary and %d '
                    'off-vocabulary value(s) emitted%s.'
                    % (len(filled), len(self.columns), n_abst, n_model, self.engine.name,
                       len(violations),
                       '; ' + ', '.join(missing) if missing else ''))
        how = ' '.join(sent)

        return dict(self._summary_fields(rows, side, filled, violations),
                    how_decided=how)

    def _summary_fields(self, rows, side, filled, violations):
        s = self.state
        g = s.get('grammar') or {}
        cov = (g.get('scored') or {}).get('coverage')
        fams = (g.get('scored') or {}).get('family_counts') or {}
        groups = s.get('family_values') or {}
        ext = s.get('external') or []
        census = s.get('census') or {}
        in_scope_exts = sorted({a['ext'] for a in s.get('in_scope', [])})
        oos = {k: v for k, v in census.items() if k not in in_scope_exts}
        heads = s.get('headers') or {}
        gate = s.get('gate') or {}
        n_model = sum(1 for r in self._log_records()
                      if r.get('decided_by') not in ('code', 'abstain'))
        n_abst = sum(1 for r in self._log_records() if r['id'].startswith('D-ABSTAIN'))
        conv = sum(1 for x in side if x.get('conversion_required') == 'true')
        searchable = sum(1 for x in side if x.get('searchable_now') == 'true')
        return {
            'accession': self.acc,
            'title': (s.get('proxi', {}) or {}).get('title', ''),
            'assay': s.get('assay', ''),
            'how_decided': '',          # filled by the caller; kept here to fix column order
            'rows_emitted': len(rows),
            'cache_paths': len(s.get('files', [])),
            'in_scope_formats': ','.join(in_scope_exts),
            'out_of_scope': '; '.join('%s=%d' % (k, v) for k, v in sorted(oos.items())),
            'grammar_coverage': '' if cov is None else round(cov, 4),
            'grammar_families': '; '.join('%s=%d' % (k, v) for k, v in sorted(fams.items())),
            'groups_decided': '; '.join(sorted(groups)),
            'submitter_table': '; '.join('%s on %s @%.3f' % (e['table'].rsplit('/', 1)[-1], e['column'],
                                                             e['result']['match_rate']) for e in ext),
            'headers_sampled': len(heads),
            'searchable_now': searchable,
            'conversion_required': conv,
            'gnps2_ms2_counted': s.get('ms2_counted', 0),
            'redu_columns': len(self.columns),
            'columns_filled': len(filled),
            'filled_column_names': ','.join(filled),
            'abstentions': n_abst,
            'cv_violations': len(violations),
            'decisions_total': self.log.n,
            'decisions_by_model': n_model,
            'engine': self.engine.name,
            'publication_gate': ('open' if gate.get('opened') else 'declined') if gate
                                else 'not reached',
            'sources_used': '; '.join(self._sources()[0]),
            'information_missing': '; '.join(self._sources()[1]),
            'run_date': __import__('datetime').date.today().isoformat(),
        }

    def _sources(self):
        s = self.state
        heads = s.get('headers') or {}
        ext = s.get('external') or []
        gate = s.get('gate') or {}
        sources = []
        if s.get('proxi'):
            sources.append('MassIVE PROXI record')
        if s.get('files'):
            sources.append('GNPS2 file cache (%d rows)' % len(s['files']))
        if heads:
            sources.append('%d sampled mzML headers' % len(heads))
        if ext:
            sources.append('%d submitter table(s)' % len(ext))
        if (s.get('ontology') or {}).get('by_uberon') or (s.get('ontology') or {}).get('by_taxonomy'):
            sources.append('OLS4 ontology lookups')
        missing = []
        if s.get('files') and not heads:
            missing.append('no readable raw-file headers')
        if s.get('files') and not ext:
            missing.append('no submitter metadata table')
        if s.get('files') and not s.get('ms2_counted'):
            missing.append('GNPS2 has not indexed the spectra')
        if gate and not gate.get('opened'):
            missing.append('publication gate declined')
        return sources, missing

    def _log_records(self):
        if not os.path.exists(self.log.path):
            return []
        out = []
        for line in open(self.log.path, encoding='utf-8'):
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except ValueError:
                    pass
        return out

    @staticmethod
    def _write_csv(path, rows, cols):
        with open(path, 'w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c, '') for c in cols})

    @staticmethod
    def _upsert_csv(path, row, key):
        """Fleet mode: one row per dataset in a shared file, replaced on re-run."""
        existing, cols = [], list(row)
        if os.path.exists(path):
            with open(path, encoding='utf-8', newline='') as f:
                rd = csv.DictReader(f)
                cols = list(dict.fromkeys(list(rd.fieldnames or []) + cols))
                existing = [r for r in rd if r.get(key) != row[key]]
        with open(path, 'w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
            w.writeheader()
            for r in existing + [row]:
                w.writerow({c: r.get(c, '') for c in cols})

    @staticmethod
    def _write_tsv(path, rows, cols):
        with open(path, 'w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c, '') for c in cols})

    def _report(self, rows, side, violations, filled):
        s = self.state
        fam = Counter(x['grammar_family'] for x in side if x.get('grammar_family'))
        conv = Counter(x['conversion_recipe'] for x in side if x['conversion_recipe'])
        L = ['# %s - metadata recovery report\n' % self.acc,
             '**Title.** %s  ' % (s['proxi'].get('title') or ''),
             '**Assay.** %s | **Engine.** %s | **Decisions logged.** %d\n'
             % (s.get('assay'), self.engine.name, self.log.n),
             '## Scope\n',
             '| | |\n|---|---|\n',
             '| cache rows | %d |\n' % len(s['files']),
             '| in-scope acquisitions | %d |\n' % len(s['in_scope']),
             '| format census | %s |\n' % json.dumps(s['census']),
             '| GNPS2 MS2-counted | %d |\n' % s['ms2_counted'],
             '\n## Grammar\n']
        g = s.get('grammar')
        if g:
            L.append('Coverage **%.1f%%** over %d acquisitions; families %s.\n\n'
                     % (100 * g['scored']['coverage'], g['scored']['n_total'], dict(fam)))
            L.append('```json\n%s\n```\n' % json.dumps(g['patterns'], indent=1))
        else:
            L.append('No grammar cleared the %.0f%% threshold; structural variables are empty.\n'
                     % (100 * self.cfg.grammar_coverage))
        if conv:
            L.append('\n## Conversion queue\n\n')
            for k, v in conv.items():
                L.append('- %d files -> %s\n' % (v, k))
        L.append('\n## Validation\n\n')
        L.append('- CV violations: **%d**\n' % len(violations))
        L.append('- filled ReDU columns: %d of %d\n' % (len(filled), len(self.columns)))
        L.append('- rows + out-of-scope = cache rows: %d + %d\n'
                 % (len(rows), len(s['files']) - len(rows)))
        L.append('\n## Filled columns\n\n' + ', '.join('`%s`' % c for c in filled) + '\n')
        L.append('\n## Gate\n\n```json\n%s\n```\n' % json.dumps(s.get('gate'), indent=1))
        return ''.join(L)

    # ------------------------------------------------------------------ all --
    def run(self, header_k=None, candidates=None):
        self.acquire()
        self.classify_assay()
        if self.state['assay'] == 'proteomics':
            self.log.record('D-ASSAY', 0, 'Continue with metabolomics recovery?', 'stop',
                            'This deposit is proteomics; metabolomics metadata recovery would pollute ReDU.',
                            [], 'A', ['whole run'])
            summary = self.summary_row([], [], [], [])
            self._write_csv(os.path.join(self.outdir, 'decision_summary_%s.csv' % self.acc),
                            [summary], list(summary))
            if self.summary_path:
                self._upsert_csv(self.summary_path, summary, 'accession')
            n = self.log.n
            if not self.keep_log and os.path.exists(self.log.path):
                os.remove(self.log.path)
            return {'stopped': 'proteomics', 'decisions': n,
                    'outputs': ['decision_summary_%s.csv' % self.acc]}
        self.sample_headers(header_k)
        self.induce_grammar()
        self.try_external_tables()
        self.ground_to_cv()
        self.per_family_cv()
        self.resolve_ontologies()
        self.gate()
        self.link_publication(candidates)
        out = self.emit()
        if not self.keep_log and os.path.exists(self.log.path):
            os.remove(self.log.path)
        return out
