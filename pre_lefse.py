import sys
import os
# import glob
import pandas as pd
from ww import f

project = sys.argv[1]

def if_exists(path):
	if not os.path.exists(path):
		os.makedirs(path)

# sample_info = pd.read_csv(f('{project}/sample_info.txt'), sep='\t')
with open (f('{project}/groupvs.txt')) as fg:
    groupvs_list = [i.strip() for i in fg]

# for file in glob.glob(f('{project}/03_Diversity/Beta/lefse/tmp/*.txt')):
# 	if re.search('L6', file):
file = f('{project}/03_Diversity/Beta/lefse/tmp/sample_info_L6.txt')
L6 = pd.read_csv(file, index_col=0, sep='\t')
for groupvs in groupvs_list:
	if_exists(f('{project}/03_Diversity/Beta/lefse/{groupvs}'))
	tmp_L6 = L6[L6['Group'].isin(groupvs.split('_vs_'))]
	tmp_file = f('{project}/03_Diversity/Beta/lefse/{groupvs}/L6.txt')
	tmp_L6.to_csv(tmp_file, sep='\t')
