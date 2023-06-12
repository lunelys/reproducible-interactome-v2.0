# reproducible-interactome-v2.0

This documentation describes the present tool to generate the database reproducible-interactome v2.0 (I), how to run this tool (II), some points to note (III), and then finishes with a few leads to improve the tool for a v3.0 (IV).


## I. Description of the tool

The present tool allows the generation of a protein-protein interaction database, aggregating several databases accessible from the PSICQUIC API (see [https://psicquic.github.io/](https://psicquic.github.io/) and [https://rdrr.io/bioc/PSICQUIC/](https://rdrr.io/bioc/PSICQUIC/)) into one csv file, eliminating explicit and implicit redundancies as specified by Marc Melkonian et al. in [his paper](https://pubmed.ncbi.nlm.nih.gov/35015827/). He is the author of the [reproducible-interactome](https://reproducible-interactome.genouest.org/) v1.0. The v2.0 was coded by Lunelys Runeshaw during a 2 months M1 internship, with Gwenael Rabut as supervisor.

**Please ensure that you have all the necessary files including 6 Python scripts and 1 mapping file**. These files should be placed in a single folder, regardless of the folder's name:

- main.py
- biogrid\_fetching.py
- psicquic\_fetching.py
- uniprotkb\_mapping.py
- cleaning\_data.py
- removing\_redundancies.py
- biogrid\_mi\_mapping.xlsx

There are 9 **parameters** the user can modify, all from the main.py file:

- taxids
- query
- max\_result
- format
- molecular\_interaction
- psicquic\_db\_to\_use,
- mi\_fetch\_descendants
- mi\_to\_exclude
- keep\_raw

The **3 (or 4) output files** will go in the folder where the 7 files are:

- dropped\_[query]\_[species]\_tab27.csv
- interactome\_[query]\_[species]\_[format].csv
- interactome\_[query]\_[species]\_[format]\_no\_redundancies.csv
- An optional file if keep\_raw = True, interactome\_[query]\_[species]\_[format]\_raw.csv

We will now go in detail about those 3 points: files (1), parameters (2), output files (3).


### 1. Files: what do they do?

- **main.py**: it is a wrapper file. It is also here that the parameters the user can modify are. There is no reason for the user to change anything else than the parameters, or anywhere else than in this script (apart from the BioGRID API key in the **biogrid\_fetching.py** script at line 84, see section II. How to run this tool?, subsection 3. Get your own BioGRID API key).

- **biogrid\_fetching.py**: this script is fetching BioGRID experimental evidences and adapt it from tab25 alike to a tab27 alike. Don't forget to modify line 84 with your own BioGRID API key before the first run!

- **psicquic\_fetching.py**: this script is fetching experimental evidence from the active services in the PSICQUIC registry (http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=STATUS).

- **uniprotkb\_mapping.py**: this script is fetching data from the Uniprot API to be used for the mapping geneID ids from the BioGRID service to uniprotkb ids.

- **cleaning\_data.py**: this script is cleaning the PPI data fetched previously:
    - removing rows that have no IDMs, no pubmed id, or no interaction identifiers,
    - getting the descendants of a list of IDMs specified by the user and that are excluded from the analysis (see **mi\_fetch\_descendants** parameter),
    - excluding IDMs specified by the user (see **mi\_fetch\_descendants** and **mi\_to\_exclude** parameters),
    - cleaning the protein names if it is an entrez gene/locuslink (it is the case for BioGRID data),
    - replacing the geneID ids from the bioGRID service to uniprotkb ids,
    - cleaning the protein names if it is an Uniprot (this removes the "uniprotkb:" prefix to keep only PXXXXX),
    - if the gene name is unclear/absent, we try to fetch the gene name from the geneID (this is done with the Uniprot API, the call is made from the **uniprotkb\_mapping.py** script, which uses the Uniprot API to generate a dictionary),
    - cleaning the gene names (similar what is done above for the protein names),
    - reordering columns if prot1 \> prot2. After this step prot1 \<= prot2 for all rows. This is necessary to check the redundancies between rows,
    - retrieving a list of the IDM MIs from the main file,
    - getting the ancestors of those IDM and if they are obsolete from the OLS API (see https://www.ebi.ac.uk/ols4),
    - and finally removing obsolete IDMs (in iRefIndex, some experimental evidences are annotated with 2 IDMs, one which is up to date and one which is obsolete). This step is kept in the code in case a similar problem occurs with other databases).

- **removing\_redundancies.py**: this script is eliminating explicit and implicit redundancies, adding the data from the redundant rows to the kept row, and finally creating 1 final csv file without any redundancies: **interactome\_[query]\_[species]\_[format]\_no\_redundancies.csv**

- **biogrid\_mi\_mapping.xlsx**: this file contains the mapping data from BioGRID format to tab27 format. Indeed, data fetched from the BioGRID API is not in tab27 format. It is therefore necessary to format the data fetched from BioGRID, which is done in the **biogrid\_fetching.py file**. This file has been modified from the one provided by BioGRID (see [https://wiki.thebiogrid.org/lib/exe/fetch.php/mi\_biogrid\_experiment\_map.xls](https://wiki.thebiogrid.org/lib/exe/fetch.php/mi_biogrid_experiment_map.xls)). Note that BioGRID also provides a file to match IDMs to BioGRID evidence code (see [https://wiki.thebiogrid.org/doku.php/psi-mi\_xml\_version\_2.5](https://wiki.thebiogrid.org/doku.php/psi-mi_xml_version_2.5))


### 2. Parameters: how to change the query?

- **taxids**. Type: list of strings of the taxonomy of the species you want to generate an interactome for If you want one species only: ['559292']. If you want several: ['559292', '9606']. If you want ALL species, ['\*']. Important:
  - it must be a taxid available to each queried service down the pipeline that requires a taxonomy id (Biogrid API, chosen PSICQUIC API, Uniprot API).
  - note that for all species, the Uniprot API call step might be VERY long, as there are a LOT of data to fetch

Example: `taxids = ['4932', '559292', '580240']`


- **query**. Type: string of a gene name of interest (to have all experimental evidence for a single protein) or None if you want to generate the full interactome of a species.

Example: `query = 'NAM7' | query = None`


- **max\_result**. Type: None to download everything, or integer to download a specific number of experimental evidence from each service (they will therefore add up). Important:
  - max\_result must be \< 10000 if you want to download a specific number of interactions (this limit comes from the BioGRID API)

Example: `max\_result = None`


- **format**. Type: string. Only 2 possibilities: tab27 or tab25, for retrocompatibility with the v1.0. Important:
  - if you use tab25, the protein-protein tag will still retrieve some other interactions

Example: `format = 'tab27'`


- **molecular\_interaction**. Type: string describing the type of interaction to fetch (protein-protein, protein-peptide, protein-nucleic acid, protein-complex, etc …) or None if you want all the type of interactions. Important:
  - Works only with protein-protein for the moment

Example: `molecular\_interaction = 'protein-protein'`


- **psicquic\_db\_to\_use**. Type: list of strings of the databases (services) you want PSICQUIC to fetch, or string 'all' to query all active services. Important:
  - case does not matter, but spelling does: if there is an error in the name, the service will not be recognized and will not be queried, silently.
  - iRefIndex is eliminated by default
  - if you are in tab27, BioGrid is eliminated by default (but the BioGRID data is fetched using the BioGRID API). If you are in tab25, the BioGRID data is fetched only using the PSICQUIC API, not from the BioGRID   API)

Example: `psicquic\_db\_to\_use = ['mint', 'intact']`


- **mi\_fetch\_descendants**. Type: list of strings. The string you will put here are the PSI-MI ontology terms you want to exclude from the start, for the IDM, in addition to all their descendants (recursively).

Example: `mi\_fetch\_descendants = ['MI:0063', 'MI:0362', 'MI:1088']`


- **mi\_to\_exclude**. Type: list of strings. The string you will put here are the PSI-MI ontology terms you want to exclude from the start, for the IDM.

Example: `mi\_to\_exclude = ['MI:0000', 'MI:0001', 'MI:0686', 'MI:0045']`


- **keep\_raw**. Type: Boolean. True if you want to have the optional file **interactome\_[query]\_[species]\_[format]\_raw.csv** in the end of the pipeline (see next subsection, 3. Output files: what’s inside?

Example: `keep\_raw = True`


### 3. Output files: what's inside?

If the **query** parameter is set to None, it will not appear in the filenames. Also, depending on the **taxids** parameter, [species] will have a different value in the filename: if one species only, the taxonomy of that species will appear in the filename, but if several, [species] becomes "\_MIXED\_SPECIES", and for all species it becomes "\_ALL\_SPECIES"

Example for a request with `query = None`, `taxids = ['9606']` and `format = 'tab25'`: interactome\_559292\_tab25\_no\_redundancies.csv

Example for a request with `query = 'NAM7'`, `taxids = ['4932', '559292', '580240']`, and `format = 'tab27'`: interactome\_NAM7\_MIXED\_SPECIES\_tab27\_no\_redundancies.csv

Example for a request with `query = None`, `taxids = ['\*']`, and `format = 'tab27'`: interactome\_ALL\_SPECIES\_tab27\_no\_redundancies.csv

- **dropped\_[query]\_[species]\_tab27.csv**: Each dropped row previously fetched from the APIs will be stored here, if the user wants to check if there is any unexpected result.

- **interactome\_[query]\_[species]\_[format].csv**: This is the file where all the experimental evidence will be stored after cleaning. The redundancies are not removed, the data not aggregated.

- **interactome\_[query]\_[species]\_[format]\_no\_redundancies.csv**: It is the file containing the non-redundant database, in a csv format. The rows that are kept are aggregated, meaning all the data from the rows that are thrown away are added to the kept rows.

- **An optional file if keep\_raw = True, interactome\_[query]\_[species]\_[format]\_raw.csv**: This is the file where all the raw experimental evidence are stored (it corresponds to the raw data fetched from the PSICQUIC and BioGRID APIs without any cleaning.


## II. How to run this tool?

_Note: I assume that if you are on Linux you need less guidance with the basics. Windows users can follow this section._


### 1. Make sure you have Python installed on your computer

- **Windows:** Open a command line tool such as Windows Terminal or Command Prompt. In the command line, type '`python`'. If Python is installed, you should see a message like "Python 3. x.x" followed by the Python prompt, which looks like this "\>\>\>". If it doesn't work, you can try '`py`' or '`python3`'.

Otherwise, you'll have to install it: from [the official python website](https://www.python.org/downloads/windows/), download the latest stable release (currently python 3.11.3)

- **Ubuntu:** To check if Python is installed on your system, open your terminal by pressing Ctrl + Alt + T. Type in "`python3`" and press Enter. If you see the Python prompt "\>\>\>", then you have Python installed on your Ubuntu machine.


### 2. Install all modules

Depending on your package manager, install the necessary dependencies. You can usually do that through the IDE, or via the command line. Pip is the most widespread package manager

- **Windows:**

In command line, go to the Python folder.

Upgrade beforehand pip: `py -m pip install --upgrade pip`

And finally install the modules: `py -m pip install requests pandas openpyxl datetime numpy`

Check [https://www.dataquest.io/blog/install-pip-windows/](https://www.dataquest.io/blog/install-pip-windows/%20)if you need help with pip. You can also check the official pip documentation on [https://pip.pypa.io/en/stable/getting-started/](https://pip.pypa.io/en/stable/getting-started/).


### 3. Get your own BioGRID API key

You need to get your own BioGRID API key before the first run at [https://webservice.thebiogrid.org/](https://webservice.thebiogrid.org/), and modify line 84 (`params['accesskey']`) of the **biogrid\_fetching.py** script with that key.


### 4. Run!

You can run the script in 2 ways:

- **Through an IDE:** some free options are Visual Studio Code, SublimeText, Idle, Vim... Hit run from the main.py script.

OR

- **Directly in command line:** go to the location of the folder where all the scripts are located, and depending on your OS and your Python version:

   - **Windows** : `cd "C:\Users\[username]\[Folder1]\[Folder2]"` and then `py main.py`. Consult [https://docs.python.org/3/faq/windows.html](https://docs.python.org/3/faq/windows.html) if you are unsure how to make it work.
   
   - **Ubuntu** : `python3 main.py`


## III. Some points to note

- iRefIndex is not used by default because its data is particularly messy.
- PCA BioGrid annotation is mapped to the IDM "protein complementation assay(MI:0090) ". This is a problem because this IDM is a parent of "Two Hybrid". Contact PSI-MI to see if they could create a PSI-MI that would enable to separate a "enzymatic complementation assay" from a "transcriptional complementation assay".
- If there are experimental evidences that do not contain a clear gene name, the corresponding rows are stored in a file named no\_gene\_name.csv. Those rows are kept in the main output file no\_redundancies, but to investigate: it would be interesting to fetch additional Uniprot data mapping left genes from other organisms (see section IV).


## IV. A few leads to improve the tool for a v3.0

- Make full use of the tab27 format, in particular with the handling of the experimental role.
- Add a request to the pubmed API with the pub\_id column to fetch authors of the publications to have perfect consistency in all the authors column.
- Add a request to the Uniprot API using the "prot1/2" columns to ensure perfect consistency in the "gene name" column.
- Add a request to the Uniprot API using the "species1/2" columns to fetch additional Uniprot data mapping left genes from other organisms if there are still no gene names at the end of the cleaning step.
- Add the possibility to query multiple gene names (instead of only one in the current version) and to input other identifiers (not only gene names, but also systematic names, aliases or uniprot protein identifiers.
- Some proteins are encoded by several genes (TEF1&TEF2, HHT1&HHT2, special isoforms case…). Think how to properly refer to aggregate experimental evidences with such proteins.
- Multiple rows have two BioGRID interaction identifiers. This shows that our method was not able to distinguish two distinct experimental evidences, that were differentially annotated by BioGRID. Look at these cases to find how to improve our aggregation method (for instance by taking into account the experimental role or the interaction type).
- Check if there are rows with two pubmed id (this was the case with iRefIndex). If yes it probably means that some of the fetched databases contains aggregated data and should be excluded from the query.
- In some cases, the protein id is missing (in the file with dropped experimental evidences, see lines under "Number of dropped experimental evidences that do not have a uniprotkb or entrez gene protein id". Look for ways to retrieve these missing information, perhaps by fetching additional data from [the Uniprot API](https://www.uniprot.org/help/api_queries) for cross references, to retrieve the prot1/2 names that do not have a uniprotkb id and are therefore thrown during the cleaning step.
- Fetch additional data from [the Uniprot API](https://www.uniprot.org/help/api_queries), like structural data (see [return fields](https://www.uniprot.org/help/return_fields), Structure and Family & Domains). Example of query with NAM7: [https://rest.uniprot.org/uniprotkb/search?query=P30771+AND+(organism\_id:559292)&fields=structure\_3d,ft\_strand,ft\_helix,ft\_turn](https://rest.uniprot.org/uniprotkb/search?query=P30771+AND+(organism_id:559292)&fields=structure_3d,ft_strand,ft_helix,ft_turn)
