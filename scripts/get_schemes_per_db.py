#!/usr/bin/env python3

import os
import sys
import argparse
import json
import subprocess
import urllib.request
from xml.etree import ElementTree

def get_db_file_content(fileUrl):
    response = urllib.request.urlopen(fileUrl)
    data     = response.read()      # a `bytes` object
    # content  = data.decode('utf-8')
    content  = json.loads(data)
    return content


mainUrl = "http://rest.pubmlst.org/db"
allDBs = get_db_file_content(mainUrl)
keyTerms = ['MLST', 'cgMLST']
excludedDbs = ['rmlst', 'test']

# for i in range(0, len(allDBs)):
#     organism = allDBs[i]['description']
#     dbs = allDBs[i]['databases']
#     for db in dbs:
#         dbSeqDefURL = db['href']
#         dbDesc = db['description']
#         if ('seqdef' in dbSeqDefURL) and not any(term in dbSeqDefURL for term in excludedDbs):
# #             print('{}\t{}\t{}\t{}'.format(i+1, organism, dbSeqDefURL))
#             schemes = get_db_file_content(dbSeqDefURL+"/schemes")['schemes']         
#             for scheme in schemes:
#                 sDesc = scheme['description']
#                 if 'cgMLST' in sDesc:# or sDesc == 'MLST':
#                     print('{}\t{}\t{}\t{}\t{}'.format(i+1, organism, dbDesc, sDesc, dbSeqDefURL))

print('{}\t{}\t{}\t{}\t{}'.format('id', 'Display Name', 'Normalized Name', 'Scheme', 'URL', 'Path'))

for i in range(0, len(allDBs)):
    organism = allDBs[i]['description']
    dbs = allDBs[i]['databases']
    for db in dbs:
        dbSeqDefURL = db['href']
        dbDesc = db['description']
        if ('seqdef' in dbSeqDefURL) and not any(term in dbSeqDefURL for term in excludedDbs):
            schemes = get_db_file_content(dbSeqDefURL+"/schemes")['schemes']         
            for scheme in schemes:
                sDesc = scheme['description']
                if 'cgMLST' in sDesc:# or sDesc == 'MLST':
                    print('{}\t{}\t{}\t{}\t{}'.format(i+1, organism, organism, sDesc, dbSeqDefURL, 'Path'))