# ReDUx

An agent that recovers ReDU-compliant sample metadata for public mass-spectrometry deposits.
**Python calls every tool; the model is consulted only where judgement is required**, and every
decision is logged with who made it.

See [`SPEC.md`](SPEC.md) for the project brief. Implementation notes:
[`metadata_agent/README.md`](metadata_agent/README.md).

```bash
python -m metadata_agent MSV000102356 --out runs/MSV000102356 --engine auto --summary runs/summary.csv
```

## Outputs — two CSVs

| file | contents |
|---|---|
| `redu_<ACC>.csv` | ReDU's live 57-column schema plus `file_format` and `indexed`, one row per acquisition |
| `decision_summary_<ACC>.csv` | one row per dataset: `how_decided` in 2-3 sentences, the sources used, what was missing, and the counts behind it |

`--summary shared.csv` upserts the dataset row into a shared file, so a fleet run ends with one row
per dataset. [`examples/summary.csv`](examples/summary.csv) shows three deposits.

## What is code and what is the model

| Code (deterministic, testable) | Model (reasoning) |
|---|---|
| fetching, per-host rate limiting, 429/500 backoff, caching | propose filename grammars |
| sidecar exclusion, vendor-directory collapse, cross-format dedupe | choose controlled-vocabulary terms |
| **scoring** a proposed regex over every file | adjudicate contradictions between sources |
| join testing under a normalization ladder and fixed thresholds | break an assay tie the ladder cannot |
| CV **membership validation** against ReDU | judge publication candidates |
| OLS4 ontology and taxonomy lookups | pick between OLS4 candidates when several match |

The model never writes a cell. It returns a proposal — a regex, a term, a verdict — and code decides
whether it holds; a rejected proposal is re-asked with the failure attached. Controlled terms are
chosen over the **whole** vocabulary with no coded shortlisting, because keyword ranking fails
exactly where it matters (for a human stool deposit, "resuspension in 50:50 MeOH:water" ranks
`water_ground` first).

## Tools, against the spec

| SPEC tool | status |
|---|---|
| USI validator | `core.UsiResolver` — resolves a USI and probes scans for MS/MS; `indexed` column reports GNPS2 parse state |
| Ontology lookup | `core.OLS4` — EBI OLS4 for UBERON terms and NCBITaxon lineage with TAXRANK-decoded ranks |
| MASST | not implemented |
| Repository access | MassIVE (PROXI + GNPS2 dataset cache). MetaboLights and Metabolomics Workbench need their own adapters |

## Assay gate

MassIVE is mostly proteomics, so the agent classifies first and stops on anything that is not
metabolomics — writing food or body-site annotation onto a proteomics deposit is the failure mode
worth preventing. The ladder, in order: the PROXI `DatasetType` keyword (mandatory on deposits
after roughly MSV000099000); proteomics-only file types (`.mzid`, `.pepXML`, `.fasta`); PTM lists;
and precursor charge from an MS2 scan (2+/3+ peptides vs singly-charged small molecules).

## Engines

`--engine auto` takes the first available: `anthropic` (needs `ANTHROPIC_API_KEY`), `claude-cli`
(`claude` on PATH), `file` (queues questions to `pending_decisions.jsonl`, reads `answers.jsonl` on
the next run), `abstain`. The deterministic pipeline runs to completion under all four — with no
model it simply abstains where judgement was needed, and records that it did.

## Thresholds are config

`core.Config`: grammar coverage `0.95`, join match rate `0.80`, injectivity `1.05`, per-host pacing
(ProteoSAFe 12 s, USI resolver 7 s — its HTTP 500 is a disguised 429 at 10 requests/minute).
Override with `--grammar-coverage` / `--join-match-rate`.

## Status

Verified end to end on `MSV000102356` (347 acquisitions, 100% grammar coverage, submitter plate map
discovered and joined at 99.4%, 0 off-vocabulary values), `MSV000086830` (226 acquisitions
published as `mzML+raw`, indexed), and `MSV000103633` (proteomics, stopped after 6 decisions).

Known gaps: no MetaboLights or Workbench adapter, no MASST tool, no web search in the publication
stage (pass candidates with `--candidates`), `.xlsx` submitter tables skipped, and no grammar
library shared across deposits.
