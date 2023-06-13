# Author: Lunelys RUNESHAW <lunelys.runeshaw@etudiant.univ-rennes.fr>
#
# This is a wrapper file. It is also here that the parameters the user can modify are. There is no reason for
# the user to change anything else than the parameters, or anywhere else than in this script, apart from the BioGRID API key in the biogrid_fetching.py script, line 84
# The parameters (taxids, query, max_result, format, molecular_interaction, psicquic_db_to_use, mi_fetch_descendants, mi_to_exclude, keep_raw) specify the query
# FIRST RUN:
# You need to fetch your own BioGRID API key at https://webservice.thebiogrid.org/, and modify line 84 (params['accesskey']) of the biogrid_fetching.py script with that key!
# Consult the documentation for more information

import os
import sys
import csv
import datetime
import biogrid_fetching
import psicquic_fetching
import uniprotkb_mapping
import cleaning_data
import removing_redundancies

# ========================== THE CODE SHOULD ONLY BE MODIFIED AT THIS LOCATION =======================================

taxids = ['4932', '559292', '580240']  # taxids = ['559292'] if you want only for one species or several. ['*'] if you want all species
query = None  # query = 'NAM7' if you want only for one protein. None if you want everything
# max_result must be < 10000 if you want to download a specific number of interactions
max_result = None  # None = download everything
format = 'tab27'  # if you use tab25, the protein-protein tag will still retrieve some other interactions
molecular_interaction = 'protein-protein'  # if you want everything, put None here. Tested only with prot-prot for the moment
# add here the db you want psicquic to fetch, all in lowercase. Be mindful to the orthographe
# iRefIndex eliminated by default + if you are in tab27, BioGrid is eliminated by default
psicquic_db_to_use = 'all'  # 'all' = query all active services, ['mint'] if you want only one db or several
# add here the MI IDM you want to eliminate from the beginning:
mi_fetch_descendants = ['MI:0063', 'MI:0362', 'MI:1088']  # themselves + their descendants will be automatically added to mi_to_exclude
mi_to_exclude = ['MI:0000', 'MI:0001', 'MI:0686', 'MI:0045']
# MI:0000(molecular interaction), MI:0001(interaction detection method)
# MI:0686(unspecified method), MI:0045(experimental interaction detection)
# MI:0063(interaction prediction), MI:0362(inference), MI:1088(phenotype-based detection assay)
keep_raw = False  # False if you want only the interactome file cleaned, True if you want to have an additional file raw with all the data before the cleaning

# ========================== ************************************************ =========================================


def file_handler(taxids, query, max_result, format):
    tab27_headers = ['prot1', 'prot2', 'gene1', 'gene2', 'idm', 'authors', 'pub_id', 'species1', 'species2',
                     'interaction_type', 'source_databases', 'interaction_identifiers', 'confidence_score',
                     'biological_role1', 'biological_role2', 'exp_role1', 'exp_role2', 'interactor_type1',
                     'interactor_type2', 'taxid_host', 'participant_id_method1', 'participant_id_method2',
                     'service_name', 'biogrid_experimental_system', 'biogrid_description', 'biogrid_type', 'throughput']
    tab25_headers = ['prot1', 'prot2', 'gene1', 'gene2', 'idm', 'authors', 'pub_id', 'species1', 'species2',
                     'interaction_type', 'source_databases', 'interaction_identifiers', 'confidence_score',
                     'service_name']
    if taxids == ['*']:
        taxids = 'ALL_SPECIES'  # just for the filename
    elif len(taxids) > 1:
        taxids = 'MIXED_SPECIES'  # just for the filename
    else:  # there is only one species in the list, we can use it in the filename
        taxids = taxids[0]
    if query is None and max_result is None:
        interactome_filename = 'interactome_all_' + taxids + '_' + format + '.csv'
    else:
        if query is None:
            interactome_filename = 'interactome_' + taxids + '_' + format + '.csv'
        else:
            interactome_filename = 'interactome_' + query + '_' + taxids + '_' + format + '.csv'
    if os.path.exists(interactome_filename) and os.path.isfile(interactome_filename):
        os.remove(interactome_filename)
        print(interactome_filename + ' deleted')
    else:
        print(interactome_filename + ' not found, is created')
    with open(interactome_filename, mode='w') as interactome_file:
        psicquic_writer = csv.writer(interactome_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        if format == 'tab25':
            headers = tab25_headers
        elif format == 'tab27':
            headers = tab27_headers
        else:
            sys.exit('The input format is wrong. Use "tab25" or "tab27"')
        psicquic_writer.writerow(headers)
    dropped_filename = interactome_filename.replace('interactome', 'dropped')
    with open(dropped_filename, mode='w') as dropped_file:
        dropped_writer = csv.writer(dropped_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        dropped_writer.writerow(headers)
    return interactome_filename


print("Starting pipeline: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
output_file = file_handler(taxids, query, max_result, format)
geneid_dict = {}
for taxid in taxids:
    if format == 'tab27':
        print("Starting to fetch BioGRID data for taxid " + taxid + ": " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
        biogrid_fetching.fetching(output_file, taxid, query, max_result, molecular_interaction)
    print("Starting to fetch PSICQUIC data for taxid " + taxid + ": " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    psicquic_fetching.fetching(output_file, taxid, query, max_result, format, molecular_interaction, psicquic_db_to_use)
    print("Starting to fetch mapping data for taxid " + taxid + ": " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    geneid_dict = geneid_dict | uniprotkb_mapping.mapping(taxid)  # to merge the dictionaries
print("Starting to clean data: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
mi_ancestors = cleaning_data.cleaning(output_file, format, molecular_interaction, geneid_dict, mi_fetch_descendants, mi_to_exclude, keep_raw)
print("Starting to remove redundancies: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
removing_redundancies.removing(output_file, mi_ancestors)
print("End pipeline: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
