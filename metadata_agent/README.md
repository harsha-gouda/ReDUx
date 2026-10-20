# metadata_agent

Per-deposit metadata recovery where **Python calls every tool and the model only makes judgement
calls**. Each decision is logged with who made it, so any cell can be traced to either a line of
code or a model answer that code then checked.

```bash
python -m metadata_agent MSV000102356 --out runs/MSV000102356 --engine auto
```

## The split

| Done by code (deterministic, testable) | Done by the model (reasoning) |
|---|---|
| fetch PROXI, dataset cache, `uniquemri`, file bodies | propose filename grammars |
| rate limiting, 429/500 backoff, on-disk caching | choose controlled-vocabulary terms |
| wrapper stripping, sidecar exclusion, vendor-dir collapse, dedupe | adjudicate source contradictions |
| format census and scope call | break an assay tie the ladder cannot |
| **scoring** a proposed regex over every file | judge publication candidates |
| join testing under the normalization ladder and thresholds | |
| CV **membership validation** | |
| the Stage 4 gate (three boolean conditions) | |
| emitting artifacts, abstention records, arithmetic checks | |

The model never writes a cell. It returns a *proposal* — a regex, a term, a verdict — and code
decides whether it holds. A rejected proposal is re-asked with the failure attached.

## Controlled vocabulary is reasoning, not matching

`analysis.py` deliberately has no shortlist helper. Term selection gets the **whole vocabulary**
plus an evidence bundle assembled from every source at once:

```
repository_record   title, summary, keywords, species block, instruments, publications
raw_file_headers    instrument models, polarities, ionization sources, analyzers, timestamps
file_layout         format census, subtrees, example filenames
submitter_tables    any deposit table that passed the join thresholds, with sample rows
grammar             the induced filename families
```

Keyword ranking was tried and removed: for a human stool deposit, the phrase *"resuspension in
50:50 MeOH:water"* ranked `water_ground` and `water_seawater` at the top. Code's only role is to
check that the returned term exists in the vocabulary.

## Engines

`--engine auto` resolves to the first available:

| engine | when | behaviour |
|---|---|---|
| `anthropic` | `ANTHROPIC_API_KEY` set | calls the API, strict JSON, retries on validator failure |
| `claude-cli` | `claude` on PATH | shells out to `claude -p --output-format json` |
| `file` | otherwise | queues questions to `pending_decisions.jsonl`, reads `answers.jsonl` on the next run |
| `abstain` | forced | every judgement call becomes an explicit abstention |

The deterministic pipeline runs to completion under **all four**. With no model, the run still
produces a valid table — it just abstains wherever judgement was required, and says so.

## Thresholds are config, not prompt

`core.Config`: grammar coverage `0.95`, join match rate `0.80`, injectivity `1.05`, per-host pacing
(ProteoSAFe 12 s, USI resolver 7 s), retries, scope lists. Override from the CLI with
`--grammar-coverage` / `--join-match-rate`.

## Outputs — two documents

```
redu_<ACC>.csv                  ReDU's live schema (57 cols) + 2 extra, one row per acquisition
decision_summary_<ACC>.csv      one row for the dataset: how it was decided, and on what
```

`--summary shared.csv` upserts that single row into a shared file, so a fleet run ends with one
row per dataset and a re-run replaces its own row rather than appending.

The ReDU CSV carries two non-ReDU columns, appended **after** the schema so a strict consumer can
read the first 57 unchanged:

| column | values | meaning |
|---|---|---|
| `file_format` | `mzML`, `mzXML`, `raw`, `d`, or `mzML+raw` | every representation of that acquisition; the row points at the searchable copy when one exists |
| `indexed` | `yes` / `ms1_only` / `no` | whether **GNPS2 has parsed** the file. `no` means not yet indexed, NOT that the file lacks MS/MS |

One injection published twice — vendor original plus converted copy, or under both `peak/` and
`ccms_peak/` — is **one row**, not two. `file_format` records the rest.

Optional, off by default: `--keep-log` writes `decisions_<ACC>.jsonl` (every decision tagged
`decided_by: code | <engine> | abstain`), `--extras` writes the per-file derived CSV.
`pending_decisions.jsonl` / `answers.jsonl` appear only under the file hand-off engine.

## Verified behaviour

On `MSV000102356`: 697 cache rows → **347 acquisitions** (347 `.mzML.scans` sidecars excluded,
2 CSV, 1 XML), grammar **100%** coverage over three families, the submitter plate map discovered
and joined automatically at **99.4%** with 1.0 rows per file, **0** CV violations, 59 decisions
(50 code, 8 model, 1 abstention).

On `MSV000102340`: stops after 11 decisions with `stopped: proteomics` — the assay ladder reads the
`DatasetType:Proteomics` keyword before any file work begins.

## Known gaps

- Stage 5 does no web search; pass candidates with `--candidates candidates.json` and the model
  judges them against the hard-evidence rule.
- ReDU column list is static in `pipeline.REDU_COLUMNS`; refresh it when ReDU's schema changes.
- `.xlsx` submitter tables are skipped (CSV/TSV only).
- No grammar library across deposits yet: each run induces from scratch.
