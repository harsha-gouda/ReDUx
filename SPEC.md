```markdown
# <Project MERCI>, Team <1>


Members: <Harsha Gouda>
Repo: <https://github.com/harsha-gouda/ReDUx>


## 1. Task
<The agent will extract sample information form three metabolomics data repositories containg 1.2 million sample files,by collecting information for various sources. Exisitic public information about this data include depositor name, title of the deposition, description, filenames, sample headers, publications linked the data deposition ID. The agent will reason what is realiable informatoin to extract biological, and anlytical information linked to file, and create a hamonized, controlled vocabolary information regarding the sample types>
Example: <ST003797 dataset from metabolights> -> <dataset contains information from Seeds and leaves of maize metabolome in 2,706 samples analyzed using GC-MS technique in Agilent 7890A instrument> #But ouput will be in a csv format


## 2. Why an agent?
<"workflow">


## 3. Tools and data (3-6 tools)

### External services the agent calls

- <Repository record (MassIVE PROXI)>: <Fetches the depositor-authored record for a dataset: title, summary, keywords, declared species, instrument, PI/submitter and publication status. This is the primary source of truth and the first thing read>. In: <accession, e.g. MSV000102356>. Out: <JSON record; also the DatasetType keyword that settles metabolomics vs proteomics>.
- <File listing (GNPS2 dataset cache, `filename` table)>: <Lists every file path the platform holds for a deposit, with size, collection and create time. Used for the format census and to find submitter tables inside the deposit>. In: <accession>. Out: <CSV of file paths and attributes>.
- <Spectrum summary (GNPS2 dataset cache, `uniquemri` table)>: <Per-run MS1/MS2 spectrum counts, mass analyzer, ionization and DDA/DIA classification, populated only for files the platform has parsed>. In: <accession>. Out: <CSV; drives the `indexed` column (yes / ms1_only / no)>.
- <File fetch (ProteoSAFe, HTTP Range)>: <Pulls the leading bytes of an mzML to read instrument model, serial, polarity, ionization and acquisition timestamp, and downloads submitter tables such as plate maps and instrument sequence files. Headers outrank prose for every analytical column>. In: <accession + file path>. Out: <first ~70 KB of the file; parsed header facts or a parsed CSV>.
- <USI validator (GNPS2 USI resolver)>: <Checks that a file is addressable and readable from its dataset USI, and probes several scans to see whether fragmentation spectra are present. Note zero spectrum counts in the cache mean "not yet indexed", not "no MS/MS" - the probe is what separates the two>. In: <USI identifier, e.g. mzspec:MSV000102356:peak/mzML/1001_Plate1_A1.mzML>. Out: <resolves yes/no, peak counts, precursor m/z and charge per scan>.
- <Ontology lookup (EBI OLS4)>: <Resolves a chosen label to an ontology identifier, and expands a taxon into its lineage with ranks decoded through TAXRANK. Identifiers are looked up, never recalled>. In: <label such as "feces", or an NCBI taxid such as 9606>. Out: <UBERON:0001988; NCBIRank/Species/Genus/Family/Order/Class/Phylum/Kingdom>.
- <ReDU schema and vocabulary>: <Reads ReDU's own column list from its dump header with a Range request, and the allowed term list for every controlled column. The schema is never hardcoded, so a column added upstream cannot silently go unfilled>. In: <column name>. Out: <57-column schema; the permitted terms for that column>.
- <MASST (planned, not implemented)>: <Given a USI, search GNPS public data and return the hits as a table>. In: <USI>. Out: <csv table of search results>.

### Internal tools the orchestrator calls (deterministic, no model)

- <Format census and scope>: <Normalises packaging suffixes, excludes index sidecars, collapses vendor directories to one acquisition, and dedupes the same injection published in several formats (mzML + raw) or collections. One injection becomes one row>. In: <file paths>. Out: <acquisition list with formats; per-extension census>.
- <Grammar scorer>: <Runs a proposed filename regex over every acquisition and reports coverage overall and per subtree, plus field cardinality. A proposal below 95% is rejected and re-asked with the failing examples attached>. In: <regex families proposed by the model>. Out: <coverage, per-family counts, parsed fields, accept/reject>.
- <Join tester>: <Tests a submitter table against the filenames under a normalization ladder (case, zero padding, delimiter swap, A01 vs A1) and accepts only at >=80% matched with <=1.05 table rows per file>. In: <filename stems, candidate table column>. Out: <match rate, injectivity, accept/reject>.
- <CV validator>: <Checks every value written against the allowed ReDU terms; an out-of-vocabulary term is refused rather than written>. In: <column, value>. Out: <pass/fail with reason>.
- <Assay classifier>: <Decides metabolomics vs proteomics before any file work, so that a proteomics deposit is never annotated as a sample. Ladder: DatasetType keyword, then proteomics-only file types (.mzid/.pepXML/.fasta), then PTM lists, then precursor charge from an MS2 scan>. In: <PROXI record, file census, charge probe>. Out: <metabolomics / proteomics / ambiguous + evidence>.
- <Decision engine>: <Carries the five judgement calls to a model and validates what comes back. Backends: Anthropic API, Claude CLI, file hand-off queue, or abstain>. In: <question + evidence bundle + allowed answers>. Out: <validated JSON decision, or an abstention>.

Data source: <MassIVE (implemented: PROXI + GNPS2 dataset cache). MetaboLights and Metabolomics Workbench are in scope for the task but still need their own adapters>. Access checked: <yes - all services above were called live during development>


## 4. Baseline (Milestone III, Tue 5 PM)
<The simplest version that runs end to end, e.g. one LLM call with
no tools, or a fixed script. Your agent will be compared to it.>


## 5. Evaluation (Milestone IV, Wed 11 AM)
- Test cases: <5-20, each with a known answer or a pass/fail check>
- Metric: <the number you will report>
- Runs per case: <model outputs vary, so more than one if you can>
- Comparison: <baseline vs. agent, model A vs. model B, ...>
- RL component (optional): <observation, actions, reward,
  baseline policy>


## 6. Demo (Milestone V, Wed 1 PM)
<The one example you will run live, and what the audience
should see.>


## 7. Risks, limits, roles
- Biggest risk: <...>
- Fallback if you are behind on Tuesday at noon: <smaller version>
- Limits: <max steps, max tokens, expected spend from your keys>
- Roles: <who owns which part>
