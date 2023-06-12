# Author: Lunelys RUNESHAW <lunelys.runeshaw@etudiant.univ-rennes.fr>
#
# This script is fetching BioGRID experimental evidences and adapt it from tab25 alike to a tab27 alike
# FIRST RUN:
# You need to fetch your own BioGRID API key at https://webservice.thebiogrid.org/, and modify line 84 
# (params['accesskey']) of the biogrid_fetching.py script with that key!

import requests
import json
import pandas as pd
import openpyxl


def make_call(base_url, params, total, start=0, max=10000):
    # Maximum number of results is limited to 10k. Paginate to retrieve everything
    print('Processing BioGrid data: ' + str(params['start']) + '/' + str(total))
    r = requests.get(base_url, params=params)
    interactions = r.json()
    # Create a hash of results by interaction identifier
    data = {}
    for interaction_id, interaction in interactions.items():
        data[interaction_id] = interaction
        # Add the interaction ID to the interaction record, so we can reference it easier
        data[interaction_id]['INTERACTION_ID'] = interaction_id
    df = pd.DataFrame.from_dict(data, orient='index')
    df.reset_index(drop=True, inplace=True)  # it takes the wrong index otherwise
    columns = [
        'ENTREZ_GENE_A',
        'ENTREZ_GENE_B',
        'OFFICIAL_SYMBOL_A',
        'OFFICIAL_SYMBOL_B',
        'EXPERIMENTAL_SYSTEM',
        'PUBMED_AUTHOR',
        'PUBMED_ID',
        'ORGANISM_A',
        'ORGANISM_B',
        'THROUGHPUT',
        'SOURCEDB',
        'INTERACTION_ID',
        'QUANTITATION',
    ]
    df = df[columns]
    return df


def biogrid_to_tab27(df_biogrid):
    # convert BioGRID format to tab27 format, with the mapping file (biogrid_mi_mapping.xlsx)
    source_databases = 'psi-mi:"MI:0463"(biogrid)'
    service_name = 'BioGrid'
    df_mapping = pd.read_excel('biogrid_mi_mapping.xlsx', index_col='BIOGRID_evidence_code')
    df_mapping.dropna(how='all')  # otherwise we get useless empty rows at the end
    df_biogrid = df_biogrid.merge(df_mapping, left_on='EXPERIMENTAL_SYSTEM', right_on='BIOGRID_evidence_code',
                                  how='left')
    df_biogrid.insert(len(df_biogrid.columns), "service_name", service_name, True)
    df_biogrid['SOURCEDB'] = source_databases
    col_ordered = ['ENTREZ_GENE_A', 'ENTREZ_GENE_B', 'OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B', 'idm', 'PUBMED_AUTHOR',
                   'PUBMED_ID', 'ORGANISM_A', 'ORGANISM_B', 'interaction_type', 'SOURCEDB', 'INTERACTION_ID',
                   'QUANTITATION', 'biological_role1_bait', 'biological_role2_prey', 'exp_role1_bait', 'exp_role2_prey',
                   'interactor_type1_bait', 'interactor_type2_prey', 'taxid_host', 'participant_id_method1_bait',
                   'participant_id_method2_prey', 'service_name', 'EXPERIMENTAL_SYSTEM', 'BIOGRID_description',
                   'BIOGRID_type', 'THROUGHPUT']
    df_biogrid = df_biogrid[col_ordered]
    # formatting the data to comply with the other DB
    df_biogrid['ENTREZ_GENE_A'] = df_biogrid['ENTREZ_GENE_A'].apply(
        lambda x: f"entrez gene/locuslink:{x}" if x != '-' else x)
    df_biogrid['ENTREZ_GENE_B'] = df_biogrid['ENTREZ_GENE_B'].apply(
        lambda x: f"entrez gene/locuslink:{x}" if x != '-' else x)
    df_biogrid['PUBMED_ID'] = df_biogrid['PUBMED_ID'].apply(lambda x: f"pubmed:{x}" if x != '-' else x)
    df_biogrid['ORGANISM_A'] = df_biogrid['ORGANISM_A'].apply(lambda x: f"taxid:{x}" if x != '-' else x)
    df_biogrid['ORGANISM_B'] = df_biogrid['ORGANISM_B'].apply(lambda x: f"taxid:{x}" if x != '-' else x)
    df_biogrid['taxid_host'] = df_biogrid['taxid_host'].apply(lambda x: f"taxid:{x}" if x != '-' else x)
    df_biogrid['INTERACTION_ID'] = df_biogrid['INTERACTION_ID'].apply(lambda x: f"biogrid:{x}" if x != '-' else x)
    df_biogrid['QUANTITATION'] = df_biogrid['QUANTITATION'].apply(lambda x: f"score:{x}" if x != '-' else x)
    return df_biogrid

# -----------------------------------------------------


def fetching(output_file, species, query, max_result, molecular_interaction):
    try:
        offset = 0
        if not max_result:
            max_result = 10000
        # Parameters outlined in the Wiki: https://wiki.thebiogrid.org/doku.php/biogridrest
        params = {
            'accesskey': 'yourBiogridAPIkeyHere',  # fetch your BioGRID API key at https://webservice.thebiogrid.org/
            'format': 'count',  # for the first run, to know when to stop
            'taxId': species,
            'start': offset,
            'max': max_result,
            'paginate': 'true',
        }
        if params["accesskey"] == 'yourBiogridAPIkey':
            print("You first need to fetch your own BioGRID API key at https://webservice.thebiogrid.org/, and modify line 84 (params['accesskey']) of the biogrid_fetching.py script with that key!")
        if species == '*':  # BioGRID API force us to remove the species param if we require all species
            del params["taxId"]
        base_url = 'https://webservice.thebiogrid.org/interactions/?'
        if query:
            # ['PMT2', 'NMD4'] For the v3.0, how to implement list of gene for query from main.py, we can do it like that
            params['geneList'] = '|'.join([query])  # Must be | separated
            params['searchNames'] = 'true'  # Search against official names
            # true to get any interaction involving EITHER gene, false: interactions between genes:
            params['includeInteractors'] = 'true'
            # true to get interactions between the gene_list’s first order interactors:
            params['includeInteractorInteractions'] = 'false'
        if molecular_interaction == 'protein-protein':  # else: we fetch ALL interactions (nucleic acid-prot interactions etc): not tested
            # evidence list to exclude if only protein-protein interactions are desired (to pick from the biogrid_mi_mapping.xlsx)
            evidence_list = ['Affinity Capture-RNA', 'Protein-RNA', 'Dosage Growth Defect', 'Dosage Lethality',
                             'Dosage Rescue', 'Synthetic Growth Defect', 'Synthetic Haploinsufficiency',
                             'Synthetic Lethality', 'Synthetic Rescue', 'Phenotypic Enhancement', 'Phenotypic Suppression',
                             'Positive Genetic', 'Negative Genetic']
            params['evidenceList'] = '|'.join(evidence_list)
            # false -> 'evidence_list' is evidence to exclude, if true -> is evidence to show
            params['includeEvidence'] = 'false'
        dataset = pd.DataFrame()
        total = requests.get(base_url, params=params).json()
        params['format'] = 'json'  # Return results in json format instead of count
        if total > 10000:
            for _ in range((total // params['max'])+1):
                df = make_call(base_url, params, total, offset, max_result)
                params['start'] = params['start'] + max_result
                dataset = pd.concat([dataset, df])
            print('Finished downloading ' + str(dataset.shape[0]) + ' interactions from BioGRID')
        else:
            dataset = make_call(base_url, params, total, offset, max_result)
        dropped_prot = dataset.loc[(dataset['ENTREZ_GENE_A'].str.match('-')) | (dataset['ENTREZ_GENE_B'].str.match('-'))]
        header = 'Number of dropped interactions that do not have (a) protein(s) name(s): ' + str(dropped_prot.shape[0])
        print(header)
        dropped_filename = output_file.replace('interactome', 'dropped')
        pd.Series([header]).to_csv(dropped_filename, mode='a', index=False, header=False)
        dropped_prot.to_csv(dropped_filename, mode='a', index=False, header=False)
        dataset = dataset.loc[~(dataset['ENTREZ_GENE_A'].str.match('-')) & ~(dataset['ENTREZ_GENE_B'].str.match('-'))]
        dataset = biogrid_to_tab27(dataset)
        print('Final number of interactions kept from BioGRID: ' + str(dataset.shape[0]))
        dataset.to_csv(output_file, mode='a', index=False, header=False)
    except TypeError:
        print('No data from BioGrid with this query.')
