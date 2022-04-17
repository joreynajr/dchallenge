import os
import pandas as pd
import requests
import numpy as np
import sys

in_fn = sys.argv[1]
out_fn = sys.argv[2]

with open(in_fn) as fr:
    link_ids = fr.readlines()
    link_ids = [x.strip() for x in link_ids]
    
disease_index = []
for i, link_id in enumerate(link_ids):

    print('Processing: {}'.format(i))
    
    # query the OLS API for term information 
    query = 'http://www.ebi.ac.uk/ols/api/terms?iri={}'.format(link_id)
    req = requests.get(query)
    json = req.json()
    
    terms = json['_embedded']['terms']
    terms = str(terms)
    
    # checking for disease within the term json str
    if terms.find('disease') != -1: 
        disease_index.append([link_id, 'disease'])

    # checking for disorder within the term json str
    elif terms.find('disorder') != -1:
        disease_index.append([link_id, 'disorder'])

    # defaulting to other when neither disease nor disorder is in json str
    else:
        disease_index.append([link_id, 'other'])
    
# write out each disease annotated link
with open(out_fn, 'w') as fw:
    for x in disease_index:
        fw.write('\t'.join(x) + '\n')
