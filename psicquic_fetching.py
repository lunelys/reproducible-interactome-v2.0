# Author: Lunelys RUNESHAW <lunelys.runeshaw@etudiant.univ-rennes.fr>
#
# This script is fetching experimental evidences from the active services in the psicquic registry:
# http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=STATUS

import requests
from urllib.request import urlopen
import xml.etree.ElementTree as ET
import pandas as pd
import datetime
from io import StringIO
import numpy as np


class PsicquicService:
    def __init__(self, name, rest_url):
        self.name = name
        self.rest_url = rest_url


def read_url(url):
    # fetch the content of a psicquic service
    try:
        file_handle = urlopen(url)
        content = file_handle.read()
        file_handle.close()
    except IOError:
        print('Cannot open URL ' + url)
        content = ''
    return content


def read_active_services_from_registry(tags):
    # return a list of all the active services from the psicquic registry at the moment of the request
    registry_active_url = 'https://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=ACTIVE&format=xml'
    if tags:
        # Note: tags only apply to the PSICQUIC registry
        # So for exemple with tags='protein-protein', we lose ChEMBL because it doesn't have the protein-protein tag,
        # but we still fetch nucleic acid-protein or small molecule-protein interactions if the DB have those tags too
        registry_active_url = registry_active_url + '&tags=' + tags
    content = read_url(registry_active_url)
    root = ET.fromstring(content)  # Create the XML reader
    xmlns = '{http://hupo.psi.org/psicquic/registry}'
    services = []
    for service in root.findall(xmlns + 'service'):
        name = service.find(xmlns + 'name')
        restUrl = service.find(xmlns + 'restUrl')
        service = PsicquicService(name.text, restUrl.text)
        services.append(service)
    return services


def query_psicquic(psicquic_service, output_file, species, interactor=None,
                max_results=None, format='tab25', tags='protein-protein'):
    # fetch data from psicquic
    psicquic_rest_url_query = psicquic_service.rest_url
    if interactor is None:
        psicquic_part_url = psicquic_rest_url_query + 'query/'
    else:
        psicquic_part_url = psicquic_rest_url_query + 'interactor/' + interactor + '?'
    if species == '*':
        psicquic_url = psicquic_part_url + species + '?&format=' + format
    else:
        psicquic_url = psicquic_part_url + 'species%3A' + species
    if max_results is None:
        psicquic_url = psicquic_url + '?&format=' + format
    else:
        psicquic_url = psicquic_url + '?firstResult=0&maxResults=' + str(max_results) + '&format=' + format
    print('\t\tURL: ' + psicquic_url)
    r = requests.get(psicquic_url)
    if r.text != 'Format not supported: tab27':
        total = r.headers['X-PSICQUIC-Count']
        if r.text:
            print('\t\tDownloading ' + total + ' experimental evidences, please wait...')
            df = pd.read_csv(StringIO(r.text), sep="\t", header=None)
            if format == 'tab27':
                df.drop(df.iloc[:, np.r_[2, 3, 15, 22:28, 29:40]], inplace=True, axis=1)
                df.insert(len(df.columns), "service_name", psicquic_service.name, True)
                df.insert(len(df.columns), "biogrid_experimental_system", '-', True)
                df.insert(len(df.columns), "biogrid_description", '-', True)
                df.insert(len(df.columns), "biogrid_type", '-', True)
                df.insert(len(df.columns), "throughput", '-', True)
            else:  # default to tab25 format
                df.drop(df.columns[[2, 3]], axis=1, inplace=True)
                df.insert(len(df.columns), "service_name", psicquic_service.name, True)
            df.columns = list(pd.read_csv(output_file).columns)
            df.to_csv(output_file, mode='a', index=False, header=False)
        else:
            print('\t\tNo experimental evidences found in the service')
    else:
        print('\t\tFormat not supported: tab27')


# -----------------------------------------------------


def fetching(output_file, species, query, max_result, format, molecular_interaction, psicquic_db_to_use):
    offset = 0
    services = read_active_services_from_registry(molecular_interaction)
    for service in services:
        if (format == 'tab27' and service.name == 'BioGrid') or service.name == 'iRefIndex':
            continue  # we do not fetch twice BioGrid data in that case ... and we eliminate by default iRefIndex
        if psicquic_db_to_use == 'all':
            print('Service: ' + service.name + ' ================================================================== ')
            print("Starting to retrieve data from: " + service.name + ', ' + datetime.datetime.now().strftime(
                "%d/%m/%Y, %H:%M:%S"))
            query_psicquic(service, output_file, species, query, max_result, format, molecular_interaction)
            print('\n')
        else:
            if service.name.lower() in psicquic_db_to_use:
                print(
                    'Service: ' + service.name + ' ================================================================== ')
                print("Starting to retrieve data from: " + service.name + ', ' + datetime.datetime.now().strftime(
                    "%d/%m/%Y, %H:%M:%S"))
                query_psicquic(service, output_file, species, query, max_result, format, molecular_interaction)
                print('\n')
