import os
import re
import sys
import glob
from ww import f
import pandas as pd


def lefse(type):
    ## 获取比较组
    with open(f('{path}/groupvs.txt')) as fg:
        groupvs_list = set([i.strip() for i in fg])
    # groups = set(sample_info['Group'].tolist())
    # groupvs_list = combinations(groups, 2)
    for file in glob.glob(f('{path}/04_Picrust/{type}/LEfSe/*.txt')):
        if re.search('lefse', file, re.I):
            ## 对所有组分析
            # os.system(f("awk -F'\t' '{{$2=null;$3=null;$4=null;$6=null;print $0}}' {file} > {path}/04_Picrust/LEfSe/L6.txt"))   #### 删除不必要的列
            # os.system(f("sed -i 's/\s\+/\t/g' {path}/04_Picrust/LEfSe/L6.txt"))    ### 替换连续空格为tab
            # os.system(f('lefse-format_input.py {path}/04_Picrust/LEfSe/L6.txt \
            #             {path}/04_Picrust/LEfSe/L6.in -f c -c 2 -s -1 -u 1 -o 1000000'))
            # os.system(f('run_lefse.py {path}/04_Picrust/LEfSe/L6.in \
            #             {path}/04_Picrust/LEfSe/L6.res -l 2'))
            # os.system(f('lefse-plot_res.py {path}/04_Picrust/LEfSe/L6.res \
            #             {path}/04_Picrust/LEfSe/L6.png --dpi 600'))
            # os.system(f('lefse-plot_cladogram.py {path}/04_Picrust/LEfSe/L6.res \
            #             {path}/04_Picrust/LEfSe/L6.cla.png --format png --dpi 600'))
            # os.system(f('lefse-plot_features.py -f diff --archive zip {path}/04_Picrust/LEfSe/L6.in \
            #             {path}/04_Picrust/LEfSe/L6.res {path}/04_Picrust/LEfSe/biomarkers.zip'))
            ## 对每个比较组分析

            data = pd.read_csv(file, index_col=0, sep='\t')
            filename = file.split('_')[-2]
            for groupvs in groupvs_list:
                tmp_group = f('{path}/04_Picrust/{type}/LEfSe/{groupvs}')
                os.system(f("mkdir -p {tmp_group}"))
                tmp = data[data['Group'].isin(groupvs.split('_vs_'))]
                tmp_file = f('{path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.txt')
                tmp.to_csv(tmp_file, sep='\t')
                os.system(f('lefse-format_input.py {path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.txt \
                            {path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.in -f c -c 2 -s -1 -u 1 -o 1000000 > /dev/null'))
                os.system(f('run_lefse.py {path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.in \
                            {path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.res -l 2 > /dev/null'))
                os.system(f('lefse-plot_res.py {path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.res \
                            {path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.png --dpi 600 > /dev/null'))

                ##lefse-plot_res.py L2.res L2.png --dpi 600 --width 20 --left_space 0.5 --feature_font_size 20 #修图参数
                # os.system(f('lefse-plot_cladogram.py {path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.res \
                #             {path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.cla.png --format png --dpi 600'))
                os.system(f('lefse-plot_features.py -f diff --archive zip {path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.in \
                            {path}/04_Picrust/{type}/LEfSe/{groupvs}/{filename}.res {path}/04_Picrust/{type}/LEfSe/{groupvs}/biomarkers.zip > /dev/null'))


path = sys.argv[1]
path = path.rstrip('/')

lefse('KEGG')
lefse('COG')

