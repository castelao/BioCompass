import os
from sys import argv
import re
import string
from collections import defaultdict, Counter
from itertools import combinations_with_replacement

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

from geneComparison import score_match

script, strain_name = argv
strain_id = os.path.basename(strain_name)

table1 = pd.read_pickle('%s.pkl' % strain_name)

#This first portion will create the distance matrix


table1['n'] = table1.index

N = table1.shape[0]
A = np.nan * np.ones(2*[N])
for idx1, idx2 in combinations_with_replacement(range(N), 2):
    gene1 = table1.iloc[idx1]
    gene2 = table1.iloc[idx2]
    score = score_match(gene1, gene2)
    A[idx1, idx2] = score
    A[idx2, idx1] = score


#This second portion will run dbscan to create a subclusters possibilities

def parse_db(db):
    D = defaultdict(list)
    for i,item in enumerate(db):
        D[item].append(i)
    D = {k:v for k,v in D.items() if len(v)>1}
    return D

def find_category(categories,col5):
    if 'biosynthetic' in categories:
        col5.append('biosynthetic')
    else:
        if len(categories) > 1:
            category = Counter(categories).most_common(1)[0][0]
            col5.append(category)
        else:
            col5.append(categories[0])
    
count = 0

db_arrays = []
for itn in range(1,len(A)):
    db = DBSCAN(eps=itn, min_samples=2).fit_predict(A)
    if not np.any([(db==d).all() for d in db_arrays]):
        db_arrays.append(db)

        subcluster_dict = parse_db(db)
        del(subcluster_dict[-1])

        output = {'BGC': [], 'subcluster': [], 'CDSs': [], 'loci': [],
                'category': []}
        for key, value in subcluster_dict.iteritems():
            output['BGC'].append(strain_name)
            output['subcluster'].append(string.ascii_uppercase[subcluster_dict.keys().index(key)])
            output['CDSs'].append(len(value))
            categories = []
            genes = []
            for item in value:
                categories.append(table1.category.loc[item])
                genes.append(table1.locus_tag.loc[item])
            genes = ','.join(genes)
            output['loci'].append(genes)
            find_category(categories, output['category'])

        count = count + 1
        table2 = pd.DataFrame(frames, index=None)
        table2.to_pickle('%s_table2_%d.pkl' % (strain_id, count))

        table2_df = pd.DataFrame(frames, index=None)
        table2_df.to_csv('%s_table2_%d.csv' % (strain_name,count), sep='\t', index=False)
