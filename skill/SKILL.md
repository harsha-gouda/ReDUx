---
name: massive-metadata-recovery
description: "Recover ReDU-compliant sample metadata for a public MassIVE/GNPS metabolomics dataset from its title, description, filenames, directory structure and raw-file headers, treating publication search as a gated fallback, restricted to the acquisition formats GNPS can actually use (mzML, mzXML, and the convertible vendor formats Thermo/Waters .raw and Agilent/Bruker .d), and emit an auditable decision log explaining every filled and every deliberately empty cell. Use when given an MSV accession and asked to extract, recover, or fill in sample metadata, or to produce a ReDU TSV."
---

# MassIVE metadata recovery

Given a MassIVE accession (MSVxxxxxxxxx), produce a ReDU-compliant per-file metadata TSV, a machine-readable evidence record, a decision log, and a human review report.

## Governing principle

The repository record — title, description, submitter/PI, directory structure, filenames, raw-file headers — is the primary source of truth. Exhaust it completely before considering any literature search. **The default for any uncertain cell is empty.** A confidently wrong value propagates into every downstream repository-scale reanalysis; a blank does not. Report precision, never fill rate.

Publications are a gated fallback (Stages 4–6). A run that never reaches Stage 5 is a normal, successful run.

**Every decision is logged — including every decision to abstain.** Because the default is an empty cell, an unexplained blank is indistinguishable from a blank the agent never considered. The decision log is what makes the output auditable and therefore usable by a curator; it is a first-class deliverable, not a debug artifact. Log at every point marked **[LOG]** below, at the moment the decision is made.

## Stage -1 — Is this deposit metabolomics at all?

MassIVE is mostly proteomics. Writing a sample type, extraction solvent or body site onto a
proteomics deposit pollutes ReDU in a way that is hard to find later, so **classify before touching
any file**. Measured against known deposits of both kinds:

1. **PROXI `DatasetType` keyword** — decisive when present. Sampling five existing deposits per band
   found it on 0/5 at MSV000070000, 078000, 085000 and 092000, and 5/5 at 099000 and 102300, so it
   became mandatory somewhere between those last two. Above that line this is a one-call answer.
2. **Proteomics-only file types** — `.mzid`, `.pepXML`, `.fasta`, `.msf`, `.dat`, `.prot.xml`. A
   single `.fasta` settles it: a sequence database has no place in a metabolomics deposit.
3. **PTM list in PROXI `modifications`** — `Carbamidomethyl`, `Oxidation`, `iTRAQ8plex` and friends.
   Supporting evidence only: a confirmed proteomics deposit (MSV000079088) still reports
   "No PTMs are included".
4. **Precursor charge from an MS2 scan** — tryptic peptides come back 2+/3+ at m/z 400-1200;
   small molecules are singly charged, usually under m/z 400. This is the rung that works on old
   deposits with no keyword and no search files.

**[LOG: D-ASSAY]** the rung that decided it and its evidence. Anything not metabolomics stops here
and emits a dataset-level record saying why — a correct, complete, zero-row run.

## File formats in scope

**Only four formats produce rows: `.mzML`, `.mzXML`, Thermo/Waters `.raw` and Agilent/Bruker `.d`.** Everything else in the deposit is counted and reported, never rowed.

| Format | In scope | Why |
|---|---|---|
| `.mzML`, `.mzXML` | **yes — searchable** | the open formats GNPS/GNPS2 analyze and index; the only two that ever carry spectrum counts in the GNPS2 cache |
| `.raw` (Thermo, Waters) | **yes — conversion required** | not readable by GNPS as deposited, but convertible to mzML with MSConvert (Waters MSe/SONAR/HDMSe/HD-DDA need Symphony) |
| `.d` (Agilent, Bruker) | **yes — conversion required** | a directory, not a file; convertible with MSConvert for Agilent MassHunter and Bruker timsTOF (`analysis.tdf`), or CompassXport for other Bruker acquisitions |
| `.wiff` (SCIEX), `.lcd`, `.mzdata.xml` | no | convertible in principle, but deliberately out of scope for this skill |
| `.cdf`, `.qgd`, `.cef`, `.h5`, `.imzML` | no | processed, GC-native or imaging output; no MS/MS a MASST-style search can use |
| `.mgf` | no | GNPS2 accepts mgf for analysis, but it is not an acquisition and is not tracked as a run |

Measured against the live GNPS2 dataset cache on 2026-10-19: of 1,602,998 cached mzML runs, 612,463 were parsed and 363,073 carry MS/MS; mzXML, 744,406 / 185,673 / 115,353. Of 1,372,519 `.raw` and 71,592 `.d` runs, **zero** have ever been parsed — which is exactly why both are rowed as conversion candidates rather than as searchable data. Being an open format is necessary but not sufficient — roughly a third of cached mzML has never been indexed.

Rules that follow from this:

- **Normalise wrappers before testing the extension, and preserve case.** `x.raw.zip` and
  `y.mzML.gz` are an acquisition plus packaging: strip `.zip`/`.gz`/`.bz2`/`.tar` repeatedly before
  deciding the format, since a single `Path.stem` leaves `x.raw` behind and mis-keys the row.
  **Lowercasing while you strip is a silent disaster**: filename grammars are case-sensitive, so
  `QC_mix_1` becomes `qc_mix_1` and every proposed regex fails. Match suffixes case-insensitively;
  return the original string.
- **A sidecar is not packaging.** `.mzML.scans`, `.wiff.scan`, `.idx` describe an acquisition but
  hold no spectra. Stripping `.scans` the way you strip `.zip` turns each sidecar into a phantom
  acquisition: in MSV000102356 that produced 694 rows for 347 injections, and the inflation looks
  entirely plausible in a summary table. Count sidecars in the census; never row them.
- **One injection is one row, across formats and collections.** A deposit that publishes the vendor
  original and its converted copy (`IS_1.raw` and `IS_1.mzML`), or the same file under both `peak/`
  and `ccms_peak/`, has not doubled its sample count. Dedupe on the stem, make the searchable copy
  the representative, and record the rest — MSV000086830 is 226 acquisitions published as
  `mzML+raw`, not 452.
- **Collapse vendor directories to one acquisition.** Waters `.raw` and Agilent `.d` are directories; the cache lists every `_FUNC001.DAT` or `AcqData/*.bin` inside them. Counting cache rows as acquisitions inflates the total several-fold.
- **Emit `.raw` and `.d` rows flagged, never as searchable.** Each carries `conversion_required = true` and `searchable_now = false`. The agent does not convert anything; it produces the work queue.
- **Record which converter each vendor row needs.** Thermo and Waters `.raw`, and Agilent `.d`, go to MSConvert; Bruker `.d` goes to CompassXport unless it is timsTOF (`analysis.tdf` present), which MSConvert handles instead. Read the vendor from the cache `instrument_vendor` field or the directory contents, and say so in the row rather than leaving the queue ambiguous.
- **Count, then abstain, on everything else.** Out-of-scope files get one dataset-level `D-FORMAT` record per format with its file count — never silent omission. A deposit that is entirely out of scope is a valid, complete run that emits zero rows and says so.

## Stage 0 — Acquire (deterministic, no inference)

Fetch and cache:

1. `https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/{MSV}` → title, description, PIs, submitter + email, species, instrument, publication field.
2. `https://datasetcache.gnps2.org/datasette/database/filename.csv?_size=max&dataset__exact={MSV}` → every file path, size, `instrument_vendor`, `instrument_model`, collection, `create_time`. **Take the full listing first, then census it by format, then filter to the four in-scope formats.** The census is a deliverable; the filtered set is what gets rows.
3. `https://datasetcache.gnps2.org/datasette/database/uniquemri.csv?_size=max&dataset__exact={MSV}` → per-run `spectra_ms1`, `spectra_ms2`, `MassAnalyzer`, `Ionization`, `classification` (DDA / DIA OR TARGETED / NO_MS2). Populated only for files GNPS2 has parsed, and only ever for mzML/mzXML. `spectra_ms2 > 0` is the one hard signal that a file holds fragmentation spectra; zero counts mean *not yet indexed*, which is not the same as *no MS/MS*.
4. Raw-file headers: fetch leading bytes of representative `.mzML`/`.mzXML` files, **sampled per directory subtree**, not all N → instrument serial, acquisition datetime, polarity, scan settings, often the LC method name.
5. Check ReDU for an existing annotation of this accession. If present, diff rather than duplicate.
6. **Read ReDU's column list from ReDU**, with a Range request on the dump header
   (`bytes=0-8192` of `https://redu.gnps2.org/dump`) — 57 columns for about 6 KB, instead of a
   ~930 MB download that times out. A hardcoded list silently drops whatever was added upstream: an
   earlier 42-column constant quietly lost the whole `NCBISuperkingdom..NCBISpecies` lineage block,
   `classification` and `MS2spectra_count`, and nothing in the run reported them missing because
   they were not in the schema to begin with.

**Service manners, learned the hard way.** ProteoSAFe rate-limits after roughly three requests and
answers 429; `massive-ftp.ucsd.edu` was unreachable throughout, so file bodies have to come through
`ProteoSAFe/DownloadResultFile` with Range headers and ~12 s spacing. The USI resolver returns an
HTTP **500 page whose body says 429** above 10 requests/minute — a burst makes healthy USIs look
broken, so pace it and treat that 500 as a retry, not a verdict.

**[LOG: D-ACQ]** which sources returned data and which failed, how many headers were sampled from which subtrees, and whether a prior ReDU annotation exists. **[LOG: D-FORMAT]** the format census: every extension present, its acquisition count after wrapper normalisation and vendor-directory collapse, and in-scope or out-of-scope with the reason. A source that failed is a constraint on everything downstream and must be visible in the log.

Instrument, polarity, acquisition mode and run date are ground truth in the headers and cache. **Never infer them from prose, and never from a publication, when a header can be read.** If headers disagree across subtrees, that is a finding — mixed-polarity and mixed-instrument datasets are common and routinely mis-annotated.

Cheap derived facts here: retention-time span → gradient length; MS1/MS2 counts → DDA vs DIA vs MS1-only; `create_time` clustering → acquisition batches.

An mzML/mzXML file whose `uniquemri` row shows `spectra_ms1 > 0` and `spectra_ms2 = 0` is MS1-only and cannot serve a fragmentation search. Record it as `searchable_now = false` with reason `ms1_only`; do not treat it as a conversion candidate, because there is nothing to convert.

**Zero counts are not evidence of absence.** A deposit GNPS2 has never parsed shows
`spectra_ms1 = spectra_ms2 = 0` for every file, which says nothing about what the files contain.
Settle it with the resolver: fetch a few scans and look for `precursor_mz > 0`, which only an MS2
scan carries. MSV000102356 reads zero across all 347 files yet resolves and returns MS/MS at scan
400 (precursor m/z 168.029) and scan 1200 (m/z 116.071). Report the two states separately —
`indexed` is about the platform, `MS/MS present` is about the data. Note also that the path form
matters: `mzspec:ACC:peak/mzML/x.mzML` resolved where the `ccms_peak/` variant did not.

## Stage 1 — Filename and directory grammar induction

This is the centerpiece. Do not parse filenames one at a time; induce a grammar and validate it by execution.

**Induce over the in-scope acquisitions only**, after wrapper normalisation. Out-of-scope files distort token cardinality and coverage — a deposit publishing both `NEG/*.raw` and `DERIVED/*.cdf` will appear to have two grammars when it has one. Coverage is reported against the in-scope denominator, and the denominator is stated.

1. **Directory first.** The path prefix is usually human-curated and often more informative than the basename (e.g. `ccms_peak/2017.AmericanGut3K.mzXMLfiles/Blanks_LCMSsolvent/` gives year, cohort and QC role). Build a subtree tree with file counts; treat each subtree as a stratum that may carry its own metadata.
2. **Profile, then propose.** Tokenize on delimiters. For each token position compute cardinality, alphabet, monotonicity, and covariance with subtree. Give the model the *statistical profile* plus ~30 stratified examples — never a raw dump of all filenames — and ask for a named grammar, e.g. `{sample_code}_{plate_well}_{injection_order}_{run_id}.mzML`.
3. **Force a regex, then score it.** Run the proposed regex over all N files and report coverage. Coverage below ~95% means the grammar is wrong or there are several grammars (very common — often one per subtree). Iterate until coverage clears, or split by subtree. **[LOG: D-GRAM]** every grammar hypothesis tried, its regex, its overall and per-subtree coverage, and why the winner beat the runners-up.
4. **Classify each induced field**: `join_key_candidate` (high cardinality, ~1:1 with files), `factor` (low cardinality — the group/condition variable), `technical` (injection order, replicate index), `noise`. **[LOG: D-FIELD]** the role assigned to each field and the statistic that drove it.
5. **Surface ambiguity, never resolve it silently.** A token like `RA1` fits both a plate well (A–H × 1–12) and a rack position. **[LOG: D-FIELD]** the competing interpretations, the evidence for each, and explicitly what new evidence would discriminate them. Leave the cell empty if nothing does.

Maintain a growing library of induced grammars keyed by submitter — lab conventions repeat heavily across one depositor's submissions, so this step becomes largely lookup over time. When a stored grammar is reused, log the reuse and re-verify its coverage on this dataset rather than assuming it still holds.

## Stage 2 — Ground title/description to controlled vocabulary

Constrained extraction, not free generation. Load ReDU's allowed terms per column and select from that list or emit `not_applicable` / `ML_import`.

- **Do not pre-filter the vocabulary by keyword.** Choosing a term is reasoning over several
  sources at once, and keyword overlap fails exactly where it matters: for a human stool deposit the
  phrase "resuspension in 50:50 MeOH:water" ranks `water_ground` and `water_seawater` top. Give the
  whole allowed list plus all the evidence, then check membership.
- **Ask per material, not per deposit, and let the submitter's own table draw the line.** One answer
  for a whole deposit is how blanks inherit study-sample metadata. Group by filename-grammar family
  first; then, where a joined table splits a family into a handful of roles, ask per role. In
  MSV000102356 the `blank` family looked uniform until the Xcalibur `Sample Type` column split it
  into 15 typed `Blank` and 16 typed `QC`, which belong to different ReDU terms.
- **Ontology identifiers are looked up, not recalled.** Resolve the chosen label through EBI OLS4
  (`/api/search` with `exact=true`, then `/ontologies/ncbitaxon/terms/<iri>/hierarchicalAncestors`
  for lineage, decoding ranks through TAXRANK). That returns `UBERON:0001988` for feces and the full
  `Homo sapiens / Homo / Hominidae / Primates / Mammalia / Chordata / Metazoa` block. Where OLS4
  cannot answer — `NCBIDivision` is an NCBI taxonomy-dump field, not a rank — abstain and say so.
- **Never emit an out-of-CV term.** "stool" instead of `feces` poisons downstream reanalysis worse than a blank does.
- Emit at **dataset scope** by default. A value propagates to file rows only when Stage 1 gives a reason it should not — e.g. a blanks subtree must not inherit `NCBITaxonomy=9606`. **[LOG: D-SCOPE]** every propagation decision and every subtree excluded from one.
- Mine the description aggressively but **literally**. Solvents, percentages, matrix words, instrument words, cohort names, n and year are frequently stated verbatim. Extract what is stated; do not elaborate. "Extracted in 90% ethanol" gives `SampleExtractionMethod=ethanol` plus a dataset-level note of the concentration — it is not licence to infer the chromatography.
- **[LOG: D-CV]** for each mapping: the source phrase quoted verbatim, the CV term chosen, and any near-miss terms rejected. Log the rejections too — a phrase with no acceptable CV term is a gap in the vocabulary worth reporting upstream.

## Stage 3 — Reconcile

Cross-check description vs headers vs filenames. Agreement raises confidence. Disagreement is a finding, not an error to smooth over: record `contradiction_flags` and report them rather than silently favouring one source. **[LOG: D-CONFLICT]** each contradiction, the value each source asserts, which source won (headers outrank prose), and whether the cell was filled or left empty as a result.

## Stage 4 — Gap check: earn the right to search

List the still-empty columns. **Proceed to Stage 5 only if ALL of these hold:**

- The missing columns are high-value subject-level biological covariates (age, sex, health status, timepoint, intervention arm) — not analytical ones, which headers already settled.
- Stage 1 produced a `join_key_candidate` with plausible per-sample cardinality. Without a candidate key, per-file enrichment is impossible, so there is nothing to look for.
- The dataset has a specific searchable handle: a distinctive title, a named cohort, or a populated `publication` field in the PROXI record.

If any condition fails, **stop and emit**. This gate is expected to decline most datasets.

**[LOG: D-GATE]** — always, whichever way it goes. Record each of the three conditions as pass/fail with its evidence, and the resulting decision. A declined gate is the single most important entry in the log, because it is what explains why a whole class of columns is empty.

## Stage 5 — Link to a publication (gated, hard evidence only)

One pass of the query ladder, pooling candidates: the accession as a literal string; the dataset title as a literal; submitter + PI + matrix + year; the depositor's publication list filtered by modality.

**Accept a link only on hard evidence:**

- the accession string appears in the paper's full text or data-availability statement (strongest, and in practice the only tier worth acting on alone); or
- the dataset title or a GNPS task ID appears in the paper; or
- methods match on ≥3 independent axes (instrument model, matrix, extraction solvent, n, year) **and** an author overlaps the depositor list.

**Author overlap alone is never sufficient** — a single prolific lab's accession co-occurs with hundreds of papers.

If nothing clears the bar after one pass: record `no_confident_link` and stop. Do not broaden the queries, lower the threshold, or pick the best of a bad set. A wrong paper is worse than no paper because it looks authoritative. Datasets legitimately have zero papers.

Distinguish `source_publication` from `reanalysis_publication` — a reanalysis paper's supplement may not describe the original samples.

**[LOG: D-LINK]** every query issued verbatim, every candidate considered with its DOI, and for each an explicit accept or reject naming the criterion that decided it. Rejected candidates matter as much as accepted ones — they show a later reader that the near-miss paper was seen and ruled out deliberately, not missed.

## Stage 6 — Join (gated, hard threshold)

Fetch the supplement and identify the per-sample table: row count within an order of magnitude of the file count, containing a high-cardinality ID-like column. Normalize it — multi-row headers, merged cells and footnote markers are the norm and are the main failure mode here.

Test each `join_key_candidate` against each ID-like column under a normalization ladder: case, zero-padding, delimiter swap, prefix/suffix strip, `A01↔A1`.

**Accept only if ≥80% of files match AND the mapping is near-injective (≤1.05 supplement rows per file).** A 40% fuzzy match is a coincidence, not a join, and is exactly how half a dataset gets silently mislabelled.

If the join fails, keep the paper's contents as dataset-level facts only and ship the TSV from Stages 0–3. Never distribute a covariate across files on a partial key.

**[LOG: D-JOIN]** every (filename field × supplement column × normalization) pair tested, with its match rate and injectivity ratio, and the accept/reject call against the thresholds. Include the near-misses — a 78% match correctly rejected is exactly the decision a reviewer needs to see.

## Stage 7 — Project, validate, report

Emit four artifacts:

1. `redu_{MSV}.tsv` — ReDU columns on the **live schema**, CV-validated, one row per in-scope
   acquisition. ReDU writes `filename` as `f.` + the path inside the deposit, with **no accession in
   it** — verified across 74,384 MassIVE rows, 100% of which take that form. Four non-ReDU columns ride alongside: `file_format` (`mzML` / `mzXML` / `raw` / `d`), `searchable_now`, `conversion_required`, and `gnps2_ms2_spectra` (from `uniquemri`, blank where unparsed).
2. `evidence_{MSV}.json` — every filled cell as `{value, source_stage, source_url_or_quote, confidence_tier, decision_ids, alternatives_considered}`. `decision_ids` cross-references the decision log, so any cell in the TSV can be traced to the reasoning that produced it.
3. `decisions_{MSV}.jsonl` — the machine-readable decision log, one record per line (format below).
4. `report_{MSV}.md` — the human-readable narrative: what was filled, what was left empty and why, which hypotheses need a human call, coverage statistics, contradiction flags.

Validate before shipping: every value in the CV; filenames matching the dataset cache exactly; **every row in one of the four in-scope formats**; **every `.raw` and `.d` row flagged `conversion_required` and never `searchable_now`, and carrying a converter**; **rows + out-of-scope census = the cache's acquisition count, with the arithmetic shown**; internal consistency (a blank with a taxonomy assigned is a bug); **and every ReDU column accounted for — each is either filled with a decision id, or carries a D-ABSTAIN record.** A column that is neither is a bug in the run, not a property of the dataset.

## Decision log

One record per decision, appended **as it is made**. Never reconstruct the log at the end — a log written afterwards reads as justification for whatever was produced rather than a record of what was actually weighed, and it will quietly omit the options that were considered and dropped.

```json
{
  "id": "D-GRAM-03",
  "stage": 1,
  "question": "Which grammar explains the Blanks_LCMSsolvent subtree?",
  "options_considered": [
    {"option": "{type}_{idx}_{well}_{rep}_{run}", "coverage": 0.99, "note": "RA1 read as plate well"},
    {"option": "{type}_{idx}_{rack}_{rep}_{run}", "coverage": 0.99, "note": "RA1 read as rack position; identical fit"}
  ],
  "chosen": "{type}_{idx}_{well}_{rep}_{run}",
  "basis": "Token alphabet A-H x 1-12 matches a 96-well plate exactly; the rack reading fits equally well but is less common in this depositor's other submissions.",
  "evidence": ["regex coverage 0.99 over 412 files in subtree", "token cardinality = 96"],
  "confidence_tier": "C",
  "unresolved": "Well vs rack not fully discriminated. A plate map in the submission, or a paper methods section, would settle it.",
  "affects": ["well_position", "plate_id"]
}
```

Required fields: `id`, `stage`, `question`, `chosen`, `basis`, `confidence_tier`, `affects`. Include `options_considered` wherever more than one reading was live, and `unresolved` wherever the decision is provisional. Keep `basis` to the reason that actually decided it — one or two sentences, not a summary of everything known.

**Decision id prefixes:** `D-ACQ` (source acquisition), `D-GRAM` (grammar hypotheses), `D-FIELD` (field role and interpretation), `D-CV` (text → controlled vocabulary), `D-SCOPE` (dataset- vs file-level propagation), `D-CONFLICT` (source disagreements), `D-GATE` (Stage 4 gate), `D-LINK` (publication accept/reject), `D-JOIN` (key pairing), `D-ABSTAIN` (deliberate empty cell), `D-FORMAT` (format census and in/out-of-scope call).

### Abstention records

**Every ReDU column left empty gets a `D-ABSTAIN` record** naming which of these applies:

- `not_applicable` — the column does not apply to this dataset (e.g. `AgeInYears` on an environmental dataset).
- `no_source` — nothing in the repository record speaks to it and the gate declined a search.
- `ambiguous` — sources speak to it but do not discriminate between readings; name the readings.
- `conflicting` — sources disagree and none outranks the others; name the conflict.
- `below_threshold` — a candidate value existed but failed an explicit threshold; give the value and the number it missed.
- `out_of_scope_format` — the file exists in the deposit but is not `.mzML`, `.mzXML` or `.raw`; name the format and the count.

This is what turns "the agent left it blank" into "the agent looked, and here is exactly what was missing" — the difference between an output a curator can act on and one they have to redo from scratch.

### Narrative section of the report

Open `report_{MSV}.md` with a short prose account of the run: what the dataset turned out to be, which sources carried the weight, the two or three judgment calls that most shaped the output, whether the publication gate opened and why, and what a human should check first. Then the tables. A curator should be able to read the first half-page and know whether to trust the TSV.

## Confidence tiers

- **A** — asserted by the repository or file header (instrument model, filename, timestamps). Auto-accept.
- **B** — stated verbatim in title/description. Auto-accept, cite the phrase.
- **C** — from validated filename/directory grammar (≥95% regex coverage). Auto-accept for structural variables; present body-site or health-type inferences to a human.
- **D** — from a confirmed publication via a threshold-passing join. Auto-accept, cite table and column.
- **E** — plausible but unsupported. **Never** written to the TSV; goes in the report as an open question.

## Derived variables to emit alongside ReDU

These are computable from files and filenames alone, do not exist in ReDU, and are often the most analytically valuable output. Emit them in a sidecar table:

- `injection_order` — monotonic run counter in filename or `create_time`. Dominant driver of signal drift.
- `acquisition_batch` — subtree, date clustering, or run-counter discontinuities. The #1 confounder in repository-scale meta-analysis.
- `plate_id`, `well_position` — from the induced well token. Edge effects, carryover, plate-map errors.
- `qc_role` — `blank_extraction` / `blank_analysis` / `pooled_QC` / `standard_mix` / `study_sample`. Conflated in ReDU's `SampleType`.
- `replicate_group`, `replicate_type` — repeated sample token with differing run index.
- `run_duration` / `gradient_length` — retention-time span from the file.
- `polarity`, `ms_levels`, `acquisition_mode` (DDA/DIA/full-scan) — headers plus MS1/MS2 counts.
- `mass_accuracy_class` / `resolution_class` — derived from instrument model.
- `file_completeness` — file count vs stated n; orphans and duplicates.
- `cohort_name` — from description or directory name.
- `file_format`, `searchable_now`, `conversion_required` — scope and readiness per acquisition.
- `conversion_recipe` — for vendor rows. Thermo `.raw` and Waters full-scan/DDA `.raw`: MSConvert; Waters MSe, SONAR, HDMSe, HD-DDA: Symphony. Agilent `.d`: MSConvert, whose scan-number export for MassHunter data has known issues worth flagging on the row. Bruker `.d`: CompassXport, except timsTOF (`analysis.tdf`), which goes to MSConvert. Output mzML, 32-bit, no zlib, peak picking with the vendor algorithm for MS levels 1-2 and the PeakPicking filter placed first.
- `filename_grammar` — the induced regex itself, stored per dataset. A reusable artifact.
- `metadata_provenance_tier` per cell; `metadata_completeness_score` per dataset; `contradiction_flags`.

## Failure modes to watch

- Mixed instruments or polarities across subtrees, papered over by a single dataset-level value.
- Filenames reused across accessions.
- A grammar that fits 95% of files because one subtree dominates the count — check coverage per subtree, not just overall.
- Supplement tables whose row count matches by coincidence.
- Blanks and pooled QCs inheriting study-sample metadata.
- A decision log written at the end rather than during the run — it will read as justification, not as a trace.
- Vendor directories counted as files: one Agilent `.d` can contribute dozens of cache rows, so a deposit of 106 acquisitions reads as 3,998 "files". Collapse to the directory before rowing, and never emit an inner `AcqData/*.bin` or `*.m/*.xml` path as an acquisition.
- `.wiff.scan` index sidecars counted as separate acquisitions in the census: they are packaging for a `.wiff`, hold no spectra, and must not double the SCIEX tally.
- Waters netCDF exports appending a function number (`…_pos01.CDF` beside `…_pos.mzML`), so the out-of-scope file looks like a separate acquisition.
- Zero `spectra_ms2` read as "no MS/MS" when it means "not yet parsed by GNPS2" — the two need different flags and different follow-up.
- Treating `.mzML.scans` as packaging rather than a sidecar, so every acquisition is counted twice.
- Counting the vendor original and its converted copy as two acquisitions.
- Lowercasing filenames while stripping wrappers, which breaks every case-sensitive grammar.
- A hardcoded ReDU column list, which cannot report the columns it does not know about.
- Answering a controlled-vocabulary question once per deposit when the deposit holds blanks, QCs and
  study samples that need different answers.
- A log that records only the decisions that produced values. The abstentions are the half that tells a curator where to spend their time.

## What previous runs changed in this procedure

Each line is a rule that exists because a run went wrong, with the deposit that exposed it.

| Deposit | What it taught |
|---|---|
| `MSV000102356` (human stool, FMD) | `.mzML.scans` sidecars inflated 347 acquisitions to 694. Zero `uniquemri` counts hid files that resolve and carry MS/MS. The submitter's Xcalibur sequence files in `other/` carried acquisition order, vial position and sample role for every injection — the richest source in the deposit, and nothing advertised them. 16 files named `blank_*` were typed `QC`, which `blank_QC` resolves without discarding either reading. |
| `MSV000086830` (burger ingredients) | The deposit publishes each injection as `.raw` and converted `.mzML`; keying on (stem, format) produced 452 rows for 226 injections. |
| `MSV000103633` (PICH proteomics) | A proteomics deposit reached the pipeline at all. The assay gate now runs before any file work and stopped it in 6 decisions. |
| `MSV000079088`, `MSV000079262`, `MSV000102340` | Older proteomics deposits carry no `DatasetType` keyword, so the ladder needs the file-type and precursor-charge rungs beneath it. |
| Whole-corpus checks | Only mzML and mzXML ever carry spectrum counts: of 1.6 M cached mzML runs 612,463 are parsed; of 1.37 M `.raw` and 71.6 K `.d` runs, **zero**. Being an open format is necessary, not sufficient — roughly a third of cached mzML has never been indexed. |

Two habits survived every run and are worth keeping explicit: **the abstention is the product**, and
**every proposal gets executed before it is believed** — a grammar that looks right at 95% coverage
is wrong until the regex has been run over all N files and the unmatched examples read.
