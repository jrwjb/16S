import os
import pandas as pd
from ww import f


def generate_file(path):
    os.system(f('mkdir -p {path}/02_OTU/krona/file_list'))
    os.system(f('mkdir -p {path}/02_OTU/taxa_tree/sample_tree'))
    ## group file
    sampleinfo = pd.read_csv(f('{path}/sample_info.txt'), sep='\t')
    groupfile = sampleinfo[['SampleID', 'Group']]
    groupfile.to_csv(f('{path}/02_OTU/taxa_tree/group.txt'), index=False, header=False, sep='\t')
    ## krona file
    otu_table = pd.read_csv(f('{path}/02_OTU/otu_table_tax.txt'), index_col=0, sep='\t')
    samples = otu_table.columns.tolist()
    for i in samples[:-1]:
        tmp_list = ['taxonomy']
        tmp_list.insert(0, i)
        data = otu_table[tmp_list]
        tax_df = data['taxonomy'].str.split(';', expand=True)
        data = data.drop('taxonomy', axis=1)
        data = pd.concat([data, tax_df], axis=1)
        data.to_csv(f('{path}/02_OTU/krona/file_list/{i}.txt'), index=False, header=False, sep='\t')


