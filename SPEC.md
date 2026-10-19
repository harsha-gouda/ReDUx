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
- <USI validator>: <Checks if the samples names are availible and searchable in the datarepository from the dataset ID>. In: <USI identifier>. Out: <Indexed/Not Indexed/>.
- <Ontology lookup>: <Based on language discription, classifes samples into controlled vocabolary as described in the Ontology index>. In: <Language/Text>. Out: <EMBL ontology code and Keyword>.
- <MASST>: <Given an USI, will search the GNPS database and dowload the MASST results using the GNPS mcp>. In: <USI>. Out: <csv table of search results>.
- <tool name>: <what it does>. In: <...>. Out: <...>.

Data source: <MassIVE, MetaboLights, and Metabolomics Workbench>. Access checked: <yes/no>


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
