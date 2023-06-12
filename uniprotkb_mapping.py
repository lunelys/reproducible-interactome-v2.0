# Author: Lunelys RUNESHAW <lunelys.runeshaw@etudiant.univ-rennes.fr>
#
# This script is fetching data from the Uniprot API to be used for the mapping geneID ids from the
# BioGRID service to uniprotkb ids.

import re
import requests
from requests.adapters import HTTPAdapter, Retry

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def get_next_link(headers):
    # because we paginate the results
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    # because we paginate the results
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        response_json = response.json()['results']
        total = response.headers["x-total-results"]
        yield response_json, total
        batch_url = get_next_link(response.headers)


# -----------------------------------------------------


def mapping(species):
    # Note: we remove protein of uncertain existence (PE5, https://www.uniprot.org/help/dubious_sequences) because they have no gene name nor geneID
    # &compressed=true? to use or is used by default?
    if species != '*':
        url = 'https://rest.uniprot.org/uniprotkb/search?query=organism_id%3A' + species + '%20NOT%20existence%3A5&fields=accession,gene_primary,xref_geneid,gene_oln&size=500&format=json'
    else:
        url = 'https://rest.uniprot.org/uniprotkb/search?query=NOT%20existence%3A5&fields=accession,gene_primary,xref_geneid,gene_oln&size=500&format=json'
    progress = 0
    geneid_dict = {}
    print('Downloading of the geneID to uniprotkb mapping data...')
    for batch, total in get_batch(url):
        for object in batch:
            if object.get('genes') and object['genes'] != [{}] and len(object['uniProtKBCrossReferences']) != 0:
                geneID_id = object['uniProtKBCrossReferences'][0]['id']
                uniprotkb_id = object['primaryAccession']
                if object.get('genes')[0].get('geneName'):
                    gene_name = object['genes'][0]['geneName']['value']
                else:
                    gene_name = '-'
                if object.get('genes')[0].get('orderedLocusNames'):
                    ordered_locus_name = object['genes'][0]['orderedLocusNames'][0]['value']
                else:
                    ordered_locus_name = '-'
                geneid_dict.update({geneID_id: dict(uniprotkb_id=uniprotkb_id, gene_name=gene_name,
                                                    ordered_locus_name=ordered_locus_name)})
        progress += len(batch)
        print(f'{progress} / {total}')
    print("Mapping data downloaded from uniprotkb: complete")
    return geneid_dict
