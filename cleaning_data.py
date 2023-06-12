# Author: Lunelys RUNESHAW <lunelys.runeshaw@etudiant.univ-rennes.fr>
#
# This script is cleaning the PPI data fetched previously:
# - removing rows that have no IDMs, no pubmed id, or no interaction identifiers,
# - getting the descendants of a list of IDMs specified by the user and that are excluded from the analysis (see
# mi_fetch_descendants parameter),
# - excluding IDMs specified by the user (see mi_fetch_descendants and mi_to_exclude parameters),
# - cleaning the protein names if it is an entrez gene/locuslink (it is the case for BioGRID data),
# - replacing the geneID ids from the bioGRID service to uniprotkb ids,
# - cleaning the protein names if it is an Uniprot, (this removes the "uniprotkb:" prefix to keep only PXXXXX)
# - if the gene name is unclear/absent, we try to fetch the gene name from the geneID (this is done with the Uniprot
# API, the call is made from the uniprotkb_mapping.py script, which uses the Uniprot API to generate a dictionary),
# - cleaning the gene names (similar what is done above for the protein names),
# - reordering columns if prot1 &gt; prot2. After this step prot1 &lt;= prot2 for all rows. This is
# necessary to check the redundancies between rows,
# - retrieving a list of the IDM MIs from the main file,
# - getting the ancestors of those IDM and if they are obsolete from the OLS API (see https://www.ebi.ac.uk/ols4),
# - and finally removing obsolete IDMs (in iRefIndex, some experimental evidences are annotated with 2 IDMs,
# one which is up to date and one which is obsolete). This step is kept in the code in case a similar problem occurs
# with other databases).

import pandas as pd
import re
import datetime
import requests

# We need both "has" and "no" to handle the "one column true but the other one is false"
# For geneid_to_uniprotkb_prot:
has_uniprotkb_equivalencies = []
no_uniprotkb_equivalencies = []
# For clean_uniprotkb_name:
regex_uniprotkb = '([OPQ][0-9][A-Z0-9]|[A-NR-Z][0-9][A-Z])[A-Z0-9][A-Z0-9][0-9]([A-Z][A-Z0-9][A-Z0-9][0-9])?'
has_clear_uniprotkb_id = []
no_clear_uniprotkb_id = []
# For get_gene_name:
no_gene_name = []


def clean_pub_id(row):
    if "pubmed" in row:  # We keep a row only if it has a pubmed identification, otherwise we ditch all
        pubmeds = []  # example, pubmed:11283351
        fields = row.split('|')
        for field in fields:
            if "pubmed" in field:
                if '.' not in field:
                    pubmeds.append(field)
                    return "|".join(pubmeds)
                else:  # it is a pubmed:10.1073 type, we don't want it at all
                    row.replace('pubmed:' + field, field)
                    return row
    else:
        return row  # the row will be deleted if there is no pubmed id


def clean_idm(row, mi_obsolete):
    # Once we know which MIs are obsolete, we come back to the file to clean it from those excess terms
    value = '-'  # if there is no MIs that are okay, we will ditch the row
    fields = row.split('|')
    for field in fields:
        mi_idm = field.split('"')[1]
        if mi_idm not in mi_obsolete:
            value = field
            return value
    return value


def to_remove_contains(dropped_filename, string, col, df, text):
    # a series of 4 generic function, could be aggregated later to avoid unnecessary repetitions in the code
    dropped_prot = df.loc[df[col].str.contains(string)]
    header = 'Number of dropped experimental evidences ' + text + str(dropped_prot.shape[0])
    print(header)
    pd.Series([header]).to_csv(dropped_filename, mode='a', index=False, header=False)
    dropped_prot.to_csv(dropped_filename, mode='a', header=False, index=False)
    df = df.loc[~(df[col].str.contains(string))]
    return df


def to_keep_contains(dropped_filename, string, col1, col2, df, text):
    # a series of 4 generic function, could be aggregated later to avoid unnecessary repetitions in the code
    dropped_prot = df.loc[
        ~(df[col1].str.contains(string)) |
        ~(df[col2].str.contains(string))]
    header = 'Number of dropped experimental evidences ' + text + str(dropped_prot.shape[0])
    print(header)
    pd.Series([header]).to_csv(dropped_filename, mode='a', index=False, header=False)
    dropped_prot.to_csv(dropped_filename, mode='a', header=False, index=False)
    df = df.loc[
        (df[col1].str.contains(string)) &
        (df[col2].str.contains(string))]
    return df


def col_remove_empty_data(dropped_filename, col, df, text):
    # a series of 4 generic functions, could be aggregated later to avoid unnecessary repetitions in the code
    dropped_prot = df.loc[(df[col].str.match('-'))]
    header = 'Number of dropped experimental evidences ' + text + str(dropped_prot.shape[0])
    print(header)
    pd.Series([header]).to_csv(dropped_filename, mode='a', index=False, header=False)
    dropped_prot.to_csv(dropped_filename, mode='a', header=False, index=False)
    df = df.loc[~(df[col].str.match('-'))]
    return df


def to_remove_list(dropped_filename, list, df, text):
    # a series of 4 generic functions, could be aggregated later to avoid unnecessary repetitions in the code
    header = 'Number of dropped experimental evidences ' + text + str(len(list))
    print(header)
    pd.Series([header]).to_csv(dropped_filename, mode='a', index=False, header=False)
    if len(list) != 0:
        dropped_prot = df.loc[df.index.isin(list)]
        dropped_prot.to_csv(dropped_filename, mode='a', header=False, index=False)
        df = df.loc[~df.index.isin(list)]
    return df


def get_psicquic_query_descendants(mi_fetch_descendants, mi_to_exclude):
    # We get the IDM's MIs of all the ancestors, but only of the IDM's MIs we have in our file, to avoid excess calculs
    for mi in mi_fetch_descendants:
        api_url = 'https://www.ebi.ac.uk/ols/api/ontologies/mi/descendants?id=' + mi + '&size=500'
        response = requests.get(api_url)
        mi_infos = response.json()['_embedded']['terms']
        mi_to_exclude.append(mi)
        for mi_info in mi_infos:
            mi_to_exclude.append(mi_info['annotation']['id'][0])
    return mi_to_exclude


def get_prot_name(row, geneid_dict, total):
    # main cleaning step for the protein columns (0 and 1)
    for prot_num in [0, 1]:
        value = row[prot_num]
        # cleaning the proteins' name if it is an entrez gene/locuslink: ------------------------------------------
        if 'entrez gene/locuslink' in value:
            value = value.split(':')[1]
            for geneid in geneid_dict:  # we try mapping the geneid to a uniprot id
                if geneid == value:
                    uniprotkb_id = geneid_dict[geneid]['uniprotkb_id']
                    value = uniprotkb_id  # P1414
                    if row.name not in no_uniprotkb_equivalencies and row.name not in has_uniprotkb_equivalencies:
                        has_uniprotkb_equivalencies.append(row.name)  # to get the index of those lines
            # we didn't find an uniprot equivalency, we ditch the whole row:
            if value == row[prot_num].split(':')[1]:
                if row.name not in no_uniprotkb_equivalencies:
                    no_uniprotkb_equivalencies.append(row.name)
                    if row.name in has_uniprotkb_equivalencies:
                        has_uniprotkb_equivalencies.remove(row.name)
        # cleaning the proteins' name if it is an uniprot: ------------------------------------------
        if 'uniprotkb' in value:
            fields = value.split('|')
            for field in fields:
                if 'uniprotkb' in field:
                    parts = field.split(':')
                    value = parts[1].split('(')[0]
                    if row.name not in no_clear_uniprotkb_id and row.name not in has_clear_uniprotkb_id:
                        has_clear_uniprotkb_id.append(row.name)
        # it was a weird uniprotkb format, we ditch the whole row
        if value == row[prot_num] and not re.search(regex_uniprotkb, row[prot_num]):
            if row.name not in no_clear_uniprotkb_id:
                no_clear_uniprotkb_id.append(row.name)
                if row.name in has_clear_uniprotkb_id:
                    has_clear_uniprotkb_id.remove(row.name)
        if row.name % 100000 == 0 or row.name == int(total) - 1:
            print('Processing prot ' + str(prot_num + 1) + ': ' + str(row.name + 1) + '/' + str(total) + '...')
        row[prot_num] = value
    return row[0], row[1]


def get_gene_name(row, geneid_dict, total):
    # main cleaning step for the gene columns (2 and 3)
    for gene_num in [2, 3]:
        value = row[gene_num]
        if '(gene name)' in value:
            fields = value.split('|')
            for field in fields:
                if '(gene name)' in field:
                    parts = field.split(':')
                    value = parts[1].split('(')[0]
        # the gene name is unclear/absent, we try to fetch the gene name from the geneID to uniprotkb mapping dictionary
        if (value == row[gene_num] and row['service_name'] != 'BioGrid') or (
                row['service_name'] == 'BioGrid' and row[gene_num] == '-'):
            for dict in geneid_dict.values():
                if row[gene_num - 2] == dict['uniprotkb_id']:
                    if dict['gene_name'] != '-':
                        value = dict['gene_name']
                    else:
                        value = dict['ordered_locus_name']
        # the gene name is still unclear/absent, we try to keep only the hgnc if it is present
        if value == row[gene_num] and row['service_name'] != 'BioGrid' and 'hgnc:' in row[gene_num]:
            fields = value.split('|')
            for field in fields:
                if 'hgnc:' in field:
                    parts = field.split(':')
                    value = parts[1].split('|')[0]  # /!\ if it is alone or at the end of the line, will crash
        # the gene name is still unclear/absent, we keep the row, but it will be empty:
        if value == row[gene_num] and row['service_name'] != 'BioGrid':
            value = '-'
            no_gene_name.append(row.name)
        if row.name % 100000 == 0 or row.name == int(total) - 1:  # TODO : Remettre à 100 pour les tests
            print('Processing gene ' + str(gene_num - 1) + ': ' + str(row.name + 1) + '/' + str(total) + '...')
        row[gene_num] = value
    return row[2], row[3]


def get_mi_idm_list(row, idm_list):
    # We get a list of all the IDM's MIs in our file
    fields = row.split('|')
    for field in fields:
        mi_idm = field.split('"')[1]
        if mi_idm not in idm_list:
            idm_list.append(mi_idm)
    return row


def get_psicquic_query_ancestors(mi_list):
    # We get the IDM's MIs of all the ancestors, but only of the IDM's MIs we have in our file, to avoid excess calculs
    mi_ancestors_dict = {}
    mi_obsolete = []
    for mi in mi_list:
        api_url = 'https://www.ebi.ac.uk/ols/api/ontologies/mi/ancestors?id=' + mi + '&size=500'
        response = requests.get(api_url)
        mi_infos = response.json()['_embedded']['terms']
        ancestors_list = []
        for mi_info in mi_infos:
            if mi_info['is_obsolete'] or mi_infos[0]['description'] == [] and mi != 'MI:0000':
                mi_obsolete.append(mi)
                continue
            elif mi_info['has_children']:
                ancestors_list.append(mi_info['annotation']['id'][0])
        if mi not in mi_obsolete:
            mi_ancestors_dict.update({mi: ancestors_list})
    return mi_ancestors_dict, mi_obsolete


# -----------------------------------------------------


def cleaning(output_file, format, molecular_interaction, geneid_dict, mi_fetch_descendants, mi_to_exclude, keep_raw):
    dropped_filename = output_file.replace('interactome', 'dropped')
    df = pd.read_csv(output_file)
    if keep_raw:
        raw_file = output_file[:-4] + '_raw' + output_file[-4:]
        df.to_csv(raw_file, index=False)  # we save it before the cleaning, as a new filename if keep_raw = True
    print('Initial number of experimental evidences: ' + str(df.shape[0]))
    df['pub_id'] = df['pub_id'].apply(clean_pub_id)  # cleaning the pubmed: we keep only the pubmed id
    if molecular_interaction == 'protein-protein' and format == 'tab27':
        df = to_keep_contains(dropped_filename, 'MI:0326', 'interactor_type1', 'interactor_type2', df,
                              'fetched from psicquic that are not a protein: ')
    df = to_keep_contains(dropped_filename, 'uniprotkb|entrez gene/locuslink', 'prot1', 'prot2', df,
                          'that do not have a uniprotkb or entrez gene protein id: ')
    df = col_remove_empty_data(dropped_filename, 'idm', df, 'that do not have an idm: ')
    df = to_keep_contains(dropped_filename, 'pubmed', 'pub_id', 'pub_id', df, 'that do not have a pubmed id: ')
    df = col_remove_empty_data(dropped_filename, 'interaction_identifiers', df, 'that do not have an interaction id: ')
    idm_to_exclude = get_psicquic_query_descendants(mi_fetch_descendants, mi_to_exclude)
    df = to_remove_contains(dropped_filename, "|".join(idm_to_exclude), 'idm', df, 'that have an idm to exclude: ')
    # cleaning the species columns from all the text:
    df['species1'] = df['species1'].apply(
        lambda x: x.split('(')[0] if ('-1' not in x) and ('-2' not in x) and ('-3' not in x) else x.split('|')[0])
    df['species2'] = df['species2'].apply(
        lambda x: x.split('(')[0] if ('-1' not in x) and ('-2' not in x) and ('-3' not in x) else x.split('|')[0])
    df['taxid_host'] = df['taxid_host'].apply(lambda x: x.split('(')[0] if ('-1' not in x) and ('-2' not in x) and
                                                                           ('-3' not in x) else x.split('|')[0])
    # reassigning wrong ontology to unspecified biological roles and experimental roles:
    df.loc[df['biological_role1'] ==
           'psi-mi:"MI:0000"(unspecified)', 'biological_role1'] = 'psi-mi:"MI:0499"(unspecified role)'
    df.loc[df['biological_role2'] ==
           'psi-mi:"MI:0000"(unspecified)', 'biological_role2'] = 'psi-mi:"MI:0499"(unspecified role)'
    df.loc[df['exp_role1'] ==
           'psi-mi:"MI:0000"(unspecified)', 'exp_role1'] = 'psi-mi:"MI:0499"(unspecified role)'
    df.loc[df['exp_role2'] ==
           'psi-mi:"MI:0000"(unspecified)', 'exp_role2'] = 'psi-mi:"MI:0499"(unspecified role)'
    # cleaning the prot parts and mapping the geneid to uniprotkb id if necessary:
    print("Starting get_prot_name: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    df[['prot1', 'prot2']] = df.apply(get_prot_name, args=(geneid_dict, df.shape[0]), axis=1, result_type="expand")
    df = to_remove_list(dropped_filename, no_uniprotkb_equivalencies, df,
                        'that do not contain a uniprotkb equivalency to their entrez gene protein id: ')
    # cleaning the proteins' name if it is an uniprot
    df = to_remove_list(dropped_filename, no_clear_uniprotkb_id, df, 'that do not contain a clear uniprotkb id: ')
    # cleaning the genes' name
    print("Starting get_gene_name: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    df[['gene1', 'gene2']] = df.apply(get_gene_name, args=(geneid_dict, df.shape[0]), axis=1, result_type="expand")
    # v3.0: if there are a lot of no_gene_name, we could try to take that array and use it again in uniprotkb mapping
    if len(no_gene_name) != 0:
        print('Note: number of experimental evidences that do not contain a clear gene name: ' + str(len(no_gene_name)))
        print('Those rows are kept in the main frame, but to investigate')
        df.loc[df.index.isin(no_gene_name)].to_csv('no_gene_name.csv', mode='w', header=False, index=False)
    # to clean the source_databases that are not formatted the same:
    df['source_databases'] = df['source_databases'].apply(lambda x: x.split('(')[0] + '(' + x.split('(')[1].lower())
    # As proteins can be filled in the database in ony order, we put them all in the same order in the line:
    # prot1 = alphanumerically inferior to prot2, so that the redundancies are all took into account
    if format == 'tab27':
        df[['biological_role1', 'biological_role2']] = df[['biological_role2', 'biological_role1']] \
            .where(df['prot1'] > df['prot2'], df[['biological_role1', 'biological_role2']].values)
        df[['exp_role1', 'exp_role2']] = df[['exp_role2', 'exp_role1']] \
            .where(df['prot1'] > df['prot2'], df[['exp_role1', 'exp_role2']].values)
        df[['interactor_type1', 'interactor_type2']] = df[['interactor_type2', 'interactor_type1']] \
            .where(df['prot1'] > df['prot2'], df[['interactor_type1', 'interactor_type2']].values)
        df[['participant_id_method1', 'participant_id_method2']] = df[
            ['participant_id_method2', 'participant_id_method1']] \
            .where(df['prot1'] > df['prot2'], df[['participant_id_method1', 'participant_id_method2']].values)
    df[['gene1', 'gene2']] = df[['gene2', 'gene1']] \
        .where(df['prot1'] > df['prot2'], df[['gene1', 'gene2']].values)
    df[['species1', 'species2']] = df[['species2', 'species1']] \
        .where(df['prot1'] > df['prot2'], df[['species1', 'species2']].values)
    df[['prot1', 'prot2']] = df[['prot2', 'prot1']] \
        .where(df['prot1'] > df['prot2'], df[['prot1', 'prot2']].values)
    idm_list = []
    mi_idm_list = df['idm'].apply(get_mi_idm_list, args=(idm_list,))
    mi_ancestors, mi_obsolete = get_psicquic_query_ancestors(idm_list)
    df['idm'] = df['idm'].apply(clean_idm, args=(mi_obsolete,))
    df = col_remove_empty_data(dropped_filename, 'idm', df, 'that have only obsolete idms: ')
    df.to_csv(output_file, index=False)
    return mi_ancestors
