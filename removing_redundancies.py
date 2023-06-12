# Author: Lunelys RUNESHAW <lunelys.runeshaw@etudiant.univ-rennes.fr>
#
# This script is eliminating explicit, implicit redundancies, and creating 1 csv files summarizing the redundancies.

import pandas as pd
import datetime
import itertools
import numpy as np
external_columns = ['only_mi_idms', 'ancestors', 'count_impl']


def clean_authors(row):
    # when merging the data, we keep only one version of authors, and the cleanest one
    values = row.split('|')
    for value in values:
        if 'et al.' in value:
            return value
        elif ' ' in value and 'et al.' not in values:
            return value
        elif value != '-' and 'et al.' not in values and ' ' not in values:
            return value
        elif row == '-':
            return value
        else:
            return values


def same_merge(x):
    # to merge the data in each field but without repetitions
    if x.name not in external_columns:  # concatenate the count_expl (cond add)?
        zipped = [str(inner).split('|') for inner in list(set(x))]
        unzipped = set(list(itertools.chain.from_iterable(zipped)))
        if '-' in unzipped and unzipped != {'-'}:
            unzipped.remove('-')
        return '|'.join(unzipped)


# -----------------------------------------------------


def removing(input_file, mi_ancestors):
    df = pd.read_csv(input_file)
    no_redundancies_file = input_file[:-4] + '_no_redundancies' + input_file[-4:]
    print('Initial number of experimental evidences: ' + str(df.shape[0]))
    print("Starting to find explicit redundancies: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    count_col = df.groupby(['prot1', 'prot2', 'idm', 'pub_id']).size().reset_index(name='count')['count']
    df = df.groupby(['prot1', 'prot2', 'idm', 'pub_id']).agg(same_merge).reset_index()
    df['count_expl'] = count_col  # add a count column for the explicit redundancies
    df['count_expl'] = df['count_expl'] - 1
    print('Number of explicit redundancies: ' + str(df['count_expl'].sum()))
    print("Starting to find implicit redundancies: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    only_mi_idms = df['idm'].apply(lambda x: x.split('"')[1])
    df.insert(len(df.columns), "only_mi_idms", only_mi_idms, True)
    df.insert(len(df.columns), "ancestors", df['only_mi_idms'].map(mi_ancestors), True)
    df.insert(len(df.columns), 'impl', np.nan, True)
    grp = df.groupby(['prot1', 'prot2', 'pub_id'])
    # could be made cleaner in the next version
    for name, group in grp:
        for index, row in group.iterrows():
            for index_int, row_int in group.iterrows():
                if row['prot1'] == row_int['prot1'] and row['prot2'] == row_int['prot2'] \
                        and (row['pub_id'] in row_int['pub_id'] or row_int['pub_id'] in row['pub_id']) \
                        and row['only_mi_idms'] in row_int['ancestors']:
                    if pd.isna(row_int['impl']) and pd.isna(row_int['impl']):
                        df.at[row.name, 'impl'] = index + 1  # just to avoid the default value, 0
                        df.at[row_int.name, 'impl'] = index + 1
                    else:
                        if pd.notna(row_int['impl']):
                            df.at[row.name, 'impl'] = row_int['impl']
                        else:
                            df.at[row_int.name, 'impl'] = row['impl']
        for index, row in group.iterrows():
            for index_int, row_int in group.iterrows():
                if row['only_mi_idms'] not in row_int['ancestors'] and pd.isna(df.at[row.name, 'impl']):
                    df.at[row.name, 'impl'] = df.at[row.name, 'interaction_identifiers']
    count_col = df.groupby(['prot1', 'prot2', 'pub_id', 'impl'], dropna=False).size().reset_index(name='count')['count']
    df = df.groupby(['prot1', 'prot2', 'pub_id', 'impl'], dropna=False).agg(same_merge).reset_index()
    df['count_impl'] = count_col
    df['count_impl'] = df['count_impl'] - 1
    print('Number of implicit redundancies: ' + str(df['count_impl'].sum()))
    df['authors'] = df['authors'].apply(clean_authors)
    df.drop(['only_mi_idms', 'ancestors', 'impl'], inplace=True, axis=1)
    # df.reindex could be made cleaner in the next version
    tab27_headers = ['prot1', 'prot2', 'gene1', 'gene2', 'idm', 'authors', 'pub_id', 'species1', 'species2',
                     'interaction_type', 'source_databases', 'interaction_identifiers', 'confidence_score',
                     'biological_role1', 'biological_role2', 'exp_role1', 'exp_role2', 'interactor_type1',
                     'interactor_type2', 'taxid_host', 'participant_id_method1', 'participant_id_method2',
                     'service_name', 'biogrid_experimental_system', 'biogrid_description', 'biogrid_type',
                     'throughput', 'count_expl', 'count_impl']
    df = df.reindex(tab27_headers, axis=1)
    print('Final number of experimental evidences, without any redundancies: ' + str(df.shape[0]))
    df.to_csv(no_redundancies_file, index=False)

