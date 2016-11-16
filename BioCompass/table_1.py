import os.path
from sys import argv
from Bio import SeqIO
import pandas as pd
import re

from BioCompass import cds_from_gbk
from BioCompass import find_category_from_product
from BioCompass import get_hits

script, gb_file, strain_name, refname, cluster_number = argv

strain_id = os.path.basename(strain_name)

# Original table1_gen
table1 = cds_from_gbk(gb_file)
table1['BGC'] = strain_id

# Original category_gen
table1 = find_category_from_product(table1)
table1.sort_values('locus_tag', axis=0, ascending=True, inplace=True)
table1.drop_duplicates(subset='locus_tag', inplace=True)
table1.reset_index(drop=True, inplace=True)

# Original table_1_extender.py
short_cluster_number = re.search(r'0*([0-9]*)',cluster_number).group(1)
file_name = '../antiSMASH_input/%s/clusterblast/cluster%s.txt' % (refname,short_cluster_number)
df = get_hits(file_name)
df.rename(columns={'query_gene': 'locus_tag'}, inplace=True)
df.rename(columns={'product': 'hit_product'}, inplace=True)
#table1 = pd.merge(table1, df, on='locus_tag', how='left')
table1 = pd.merge(table1, df, on='locus_tag', how='inner')
table1['best_hit_%cov'] = table1.coverage
table1['best_hit_%id'] = table1.identity
table1['best_hit_BGC'] = table1.locus
table1['best_hit_gene'] = table1.subject_gene
table1['best_hit_gene_loc'] = table1.hit_gene


gb_record = SeqIO.read(open(gb_file, "r"), "genbank")

table1_handle = open('%s_table1.csv' % strain_name, "w")
table1.to_csv(table1_handle, sep='\t', index=False)
table1_handle.close()

gb_record.id = '%s' % strain_id

output_handle = open("%s.gbk" % strain_name, "w")
SeqIO.write(gb_record, output_handle, "genbank")
output_handle.close()