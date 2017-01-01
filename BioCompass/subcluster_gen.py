import os
from sys import argv
import string
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

def dominant_category(categories):
    if 'biosynthetic' in categories:
        return 'biosynthetic'
    return categories.mode()[0]
    
count = 0

db_arrays = []
for itn in range(1,len(A)):
    db = DBSCAN(eps=itn, min_samples=2).fit_predict(A)
    if not np.any([(db==d).all() for d in db_arrays]):
        db_arrays.append(db)

        output = {'BGC': [], 'subcluster': [], 'CDSs': [], 'loci': [],
                'category': [], 'DBSCAN_eps': []}

        # Ignores clusters -1
        for i, cluster_id in enumerate(sorted(set(db[db>=0]))):
            output['BGC'].append(strain_id)
            cluster_table = table1[db==cluster_id]
            output['CDSs'].append(len(cluster_table))
            output['loci'].append(",".join(cluster_table.locus_tag.tolist()))
            output['category'].append(
                dominant_category(cluster_table.category))
            # Is this a good idea? We're limited to 26, and with dataframe I'm not sure how subcluster would be necessary.
            output['subcluster'].append(string.ascii_uppercase[i])
            output['DBSCAN_eps'].append(itn)

        count = count + 1
        table2 = pd.DataFrame(output, index=None)
        table2.to_pickle('%s_table2_%d.pkl' % (strain_name, count))
