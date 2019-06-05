import os
import sys
import glob
import re
from ww import f

def heatmap(path):
    for file in glob.glob(path + '/02_OTU/taxa_heatmap/*/*.csv'):
        if re.search(r'group', file, re.I):
            os.system(f('Rscript /home/jbwang/code/heatmap_group.R {path}/sample_info.txt {file}'))
        else:
            os.system(f('Rscript /home/jbwang/code/heatmap.R {path}/sample_info.txt {file}'))

def venn(path):
    os.system(f('Rscript /home/jbwang/code/venn.R {path}'))

def barplot(path):
    for file in glob.glob(path + '/02_OTU/Top10*/*/*.csv'):
        os.system(f('Rscript /home/jbwang/code/top10_barplot.R {file}'))

def alpha(path):
    os.system(f('Rscript /home/jbwang/code/AlphaDiversity.R {path}'))
    for file in glob.glob(path + '/03_Diversity/Alpha/alpha_div_collated/*_T.txt'):
        os.system(f('Rscript /home/jbwang/code/rareplot.R {path} {file}'))

def beta(path):
    pcoa = f('{path}/03_Diversity/Beta/PCoA')
    tree_barplot = f('{path}/03_Diversity/Beta/Tree_barplot')
    os.system(f('mkdir {pcoa}'))
    os.system(f('mkdir -p {tree_barplot}/weighted_unifrac'))
    os.system(f('mkdir {tree_barplot}/unweighted_unifrac'))
    os.system(f('cp {path}/03_Diversity/Beta/weighted_unifrac_dm.txt {tree_barplot}/weighted_unifrac'))
    os.system(f('cp {path}/03_Diversity/Beta/unweighted_unifrac_dm.txt {tree_barplot}/unweighted_unifrac'))
    os.system(f('Rscript /home/jbwang/code/BetaDiversity.R {path}'))
    for file in glob.glob(path + '/03_Diversity/Beta/Tree_barplot/*/*_dm.txt'):
        os.system(f('Rscript /home/jbwang/code/upgma.R {path} {file}'))

def main():
    path = sys.argv[1].rstrip('/')
    heatmap(path)
    venn(path)
    barplot(path)
    alpha(path)
    beta(path)

if __name__=='__main__':
    main()
    

