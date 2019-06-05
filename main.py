import os
import re
import glob
import argparse
from itertools import combinations
from ww import f
import pandas as pd
import q30
import otu
from q30 import read_count
from q30 import stat
from raredata import redo_data
from tax import tax


class Pipeline(object):
    def __init__(self, project):
        self.project = project

    def if_exists(self, path):
        if not os.path.exists(path):
            os.makedirs(path)

    def sample_id(self):
        sample_info = os.path.join(self.project, 'sample_info.txt')
        df = pd.read_csv(sample_info, sep='\t')
        id_list = '\t'.join(list(df['SampleID']))
        return id_list

    def flash(self):
        print('开始组装...')
        self.if_exists(self.project + '/00_RawData/ExtendData')
        # for file in os.listdir(self.project + '/00_RawData'):
        # 	if re.search(r'.1.fq.gz|_1.fq.gz|.1.fastq.gz|_1.fastq.gz|.R1.fq.gz|_R1.fq.gz |.R1.fastq.gz|_R1.fastq.gz', file):
        # 		sample = '_'.join(os.path.basename(file).split('_')[:-1])
	       #      tmp = sample.split('_')[-1]
	       #      os.system(f('/home/jbwang/soft/FLASH-1.2.11/flash \
	       #                  {self.project}/00_RawData/{sample}_1.fq.gz {self.project}/00_RawData/{sample}_2.fq.gz \
	       #                  -d {self.project}/00_RawData/ExtendData -o {tmp} -x 0.1 -M 150 -t 16 -z'))

        for file in glob.glob(self.project + '/00_RawData/*_1.fq.gz'):
            sample = '_'.join(os.path.basename(file).split('_')[:-1])
            tmp = sample.split('_')[-1]
            os.system(f('/home/jbwang/soft/FLASH-1.2.11/flash \
                        {self.project}/00_RawData/{sample}_1.fq.gz {self.project}/00_RawData/{sample}_2.fq.gz \
                        -d {self.project}/00_RawData/ExtendData -o {tmp} -x 0.1 -M 150 -t 16 -z'))
        print('组装完成')

    def qiime_qc(self):
        print('开始质控及去嵌合体...')
        self.if_exists(self.project + '/00_RawData/TrimQC')
        self.if_exists(self.project + '/01_CleanData')
        self.if_exists(self.project + '/tmp')
        for file in glob.glob(self.project + '/00_RawData/ExtendData/*.extendedFrags.fastq.gz'):
            id_match = re.search(file.split('/')[-1].split('.')[0], self.sample_id()).group()
            os.system(f('split_libraries_fastq.py  \
                        -i {file} --sample_ids {id_match} \
                        -o {self.project}/00_RawData/TrimQC/{id_match} \
                        -q 19 --max_bad_run_length 3 \
                        --min_per_read_length_fraction 0.75 \
                        --max_barcode_errors 0 \
                        --store_demultiplexed_fastq \
                        --barcode_type not-barcoded'))
            os.system(f('mv {self.project}/00_RawData/TrimQC/{id_match}/seqs.fna \
                            {self.project}/00_RawData/TrimQC/{id_match}/{id_match}.fa'))
            os.system(f('mv {self.project}/00_RawData/TrimQC/{id_match}/seqs.fastq \
                            {self.project}/00_RawData/TrimQC/{id_match}/{id_match}.fq'))
            # 对每个样品去重
            os.system(f('/home/jbwang/soft/Usearch/Usearch11 -fastx_uniques {self.project}/00_RawData/TrimQC/{id_match}/{id_match}.fa \
                        -fastaout {self.project}/00_RawData/TrimQC/{id_match}/{id_match}_norep.fa'))
            # Usearch检测嵌合体
            os.system(f('/home/jbwang/soft/Usearch/Usearch11 \
                        -uchime2_ref {self.project}/00_RawData/TrimQC/{id_match}/{id_match}_norep.fa \
                        -db /home/jbwang/refedata/gold/gold.fa \
                        -chimeras {self.project}/tmp/chimeras.fa \
                        -notmatched {self.project}/tmp/notgold.fa \
                        -uchimeout {self.project}/tmp/gold.uchime \
                        -strand plus \
                        -mode balanced\
                        -threads 16'))

            # 获得嵌合体序列ID
            os.system(f("grep '>' {self.project}/tmp/chimeras.fa| sed 's/>//g' > {self.project}/tmp/chimeras.id"))
            # 去除嵌合体
            self.if_exists(self.project + f('/01_CleanData/{id_match}'))
            os.system(f('filter_fasta.py -f {self.project}/00_RawData/TrimQC/{id_match}/{id_match}.fa \
                        -o {self.project}/01_CleanData/{id_match}/{id_match}.fa \
                        -s {self.project}/tmp/chimeras.id -n'))
            os.system(f("sed -i 's/ .*//g' {self.project}/01_CleanData/{id_match}/{id_match}.fa"))
            # 获得所有序列的ID
            os.system(f("grep '>' {self.project}/00_RawData/TrimQC/{id_match}/{id_match}.fa | sed 's/>//g' > {self.project}/tmp/ids"))
            # 去除嵌合体的ID
            os.system(f("grep -F -v -f {self.project}/tmp/chimeras.id {self.project}/tmp/ids | sed 's/>//g' > {self.project}/tmp/id"))
            # 筛选fastq
            os.system(f('/home/jbwang/soft/seqtk/seqtk subseq \
                {self.project}/00_RawData/TrimQC/{id_match}/{id_match}.fq \
                {self.project}/tmp/id > {self.project}/01_CleanData/{id_match}/{id_match}.fq'))
            os.system(f('rm {self.project}/tmp/chimeras.id {self.project}/tmp/chimeras.fa {self.project}/tmp/notgold.fa \
                        {self.project}/tmp/gold.uchime {self.project}/tmp/ids {self.project}/tmp/id'))
            os.system(f('rm -rf {self.project}/01_CleanData/{id_match}/sed*'))
        print('质控完成')

    def get_fq_file(self):
        # 获取fastq文件
        rawPath = self.project + '/00_RawData'
        cleanPath = self.project + '/01_CleanData'
        rawFile_list = []
        cleanFile_list = []
        for file in os.listdir(rawPath):
            if re.search(r'_1', file):
                rawData = rawPath + '/' + file
                rawFile_list.append(rawData)
            # elif re.search(r'extend', file):
            #     combineData = data + '\\' + file
        for file in glob.glob(cleanPath + '/*/*.fq'):
            cleanFile_list.append(file)

        return sorted(rawFile_list), sorted(cleanFile_list)

    def data_statis(self):
        print('开始质控前后数据统计...')
        static_file = self.project + '/static.csv'
        with open(static_file, 'w') as f:
            f.write('SampleID,Raw PE,Clean PE,Base(nt),AvgLen(nt),Q20(%),Q30(%),GC(%),Effective(%)\n')
            for rawFile, cleanFile in zip(self.get_fq_file()[0], self.get_fq_file()[1]):
                sample = cleanFile.split('/')[-2]
                rawReads = str(read_count(rawFile))
                cleanReads = str(read_count(cleanFile))
                tmp = stat(cleanFile)
                total_count = str(tmp[0])
                average_len = str(tmp[1])
                q20_percents = str(tmp[2])
                q30_percents = str(tmp[3])
                GC_percents = str(tmp[4])
                Effect = str(round(100 * float(cleanReads) / float(rawReads), 2))

                f.write('{},{},{},{},{},{},{},{},{}\n'.format(
                    sample, rawReads, cleanReads, total_count, average_len, q20_percents, q30_percents, GC_percents, Effect))
        print('统计完成')

    def derep_fa(self):
        # 去重复
        os.system(f('cat {self.project}/01_CleanData/*/*.fna > {self.project}/01_CleanData/all.fa'))
        os.system(f("sed -i 's/_.*//g' {self.project}/01_CleanData/all.fa"))
        os.system(f('/home/jbwang/soft/Usearch/Usearch11 \
                    -fastx_uniques {self.project}/01_CleanData/all.fa \
                    -sizeout \
                    -relabel Uniq \
                    -fastaout {self.project}/01_CleanData/all_norep.fa'))

    def otu(self):
        print('开始生成OTU序列和OTU_table...')
        self.if_exists(self.project + '/02_OTU')
        # Usearch 方法
        # 聚类生成 OTU 序列
        os.system(f('/home/jbwang/soft/Usearch/Usearch11 \
                    -cluster_otus {self.project}/01_CleanData/all_norep.fa \
                    -relabel OTU_ \
                    -uparseout {self.project}/01_CleanData/uparse.txt \
                    -otus {self.project}/02_OTU/otus.fa'))
        # 生成 Otu_table
        os.system(f('/home/jbwang/soft/Usearch/Usearch11 \
                            -usearch_global {self.project}/01_CleanData/all.fa \
                            -db {self.project}/02_OTU/otus.fa \
                            -otutabout {self.project}/02_OTU/otu_table.txt \
                            -strand plus \
                            -id 0.97'))
        print('OTU完成')


    def annotation(self):
        # # 分类
        # Greengene base(默认)
        print('使用Greengene注释...')
        os.system('assign_taxonomy.py -i {}/02_OTU/otus.fa \
        -m rdp -c 0.8 -o {}/02_OTU'.format(self.project, self.project))

        # Silva base
        # print('开始注释(使用silva数据库)...')
        # os.system(f('assign_taxonomy.py -i {self.project}/OTU/otus.fa \
        #         -r /home/jbwang/refedata/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna \
        #         -t /home/jbwang/refedata/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt \
        #         -m rdp -c 0.8 -o {self.project}/OTU --rdp_max_memory 40000'))

        # otu_table.txt 转换为 biom 格式
        os.system(f('biom convert -i {self.project}/02_OTU/otu_table.txt \
                    -o {self.project}/02_OTU/otu_table.biom \
                    --table-type="OTU table" \
                    --to-json'))

        # # 添加物种信息至OTU表最后一列，命名为taxonomy
        os.system(f('biom add-metadata -i {self.project}/02_OTU/otu_table.biom \
                    --observation-metadata-fp {self.project}/02_OTU/otus_tax_assignments.txt \
                    -o {self.project}/02_OTU/otu_table_tax.biom \
                    --sc-separated taxonomy \
                    --observation-header OtuId,taxonomy'))

        os.system(f('biom convert -i {self.project}/02_OTU/otu_table_tax.biom \
                    -o {self.project}/02_OTU/otu_table_tax.txt \
                    --table-type="OTU table" \
                    --header-key taxonomy \
                    --to-tsv'))

        # # 查看OTU表的基本信息：样品，OTU数量统计
        os.system(f('biom summarize-table \
            -i {self.project}/02_OTU/otu_table_tax.biom \
            -o {self.project}/02_OTU/otu_table_tax.sum'))
        print('注释完成')

    def core_microbiome(self):
        os.system(f('compute_core_microbiome.py -i {self.project}/02_OTU/otu_table_tax.sum -o {self.project}/OTU/Core'))

    def build_tree(self):
        print('开始建树...')
        self.if_exists(self.project + '/02_OTU/backup')
        # clustalo软件比对
        os.system(f('/home/jbwang/soft/clustalo/clustalo-1.2.4-Ubuntu-x86_64 \
                    -i {self.project}/02_OTU/otus.fa \
                    -o {self.project}/02_OTU/backup/otus.afa \
                    --seqtype=DNA \
                    --full --force'))

        # muscle软件比对
        # os.system(f'/home/jbwang/soft/muscle/muscle3.8.31_i86linux64 \
        #             -in {self.project}/02_OTU/otus.fa \
        #             -out {self.project}/02_OTU/tmp/otus.afa')

        # 筛选结果中保守序列和保守区
        os.system(f('filter_alignment.py -i {self.project}/02_OTU/backup/otus.afa -o {self.project}/02_OTU/backup/'))

        # 基于fasttree建树
        os.system(f('make_phylogeny.py \
                    -i {self.project}/02_OTU/backup/otus_pfiltered.fasta \
                    -o {self.project}/02_OTU/otus.tre'))
        print('建树完成')

    def seqs_min(self):
        # 获取最小数据量
        with open(f('{self.project}/02_OTU/otu_table_tax.sum')) as f1:
            txt = f1.read()
            s = re.search('Min:.*', txt).group().split(':')[1].split('.')[0]
            min_num = int(''.join(s.split(',')))
        return min_num

    def diversity(self):
        print('开始qiime1多样性分析...')
        if os.path.exists(f('{self.project}/03_Diversity')):
            os.system(f('rm -rf {self.project}/03_Diversity'))
        # 添加#到sample_info.txt中
        os.system(f("sed -i 's/SampleID/#SampleID/' {self.project}/sample_info.txt"))
        ## -e 参数 最小数据量
        os.system(f('core_diversity_analyses.py \
                    -o {self.project}/03_Diversity \
                    -i {self.project}/02_OTU/otu_table_tax.biom \
                    -m {self.project}/sample_info.txt \
                    -t {self.project}/02_OTU/otus.tre \
                    -e {self.seqs_min()} \
                    -p /home/jbwang/code/config/alpha_params.txt'))

        # 用完修改回来
        os.system(f("sed -i 's/#//' {self.project}/sample_info.txt"))
        print('qiime1多样性分析完成')

    def alpha_diversity(self):
        # 计算常用的Alpha多样性指数
        os.system(f('mv {self.project}/03_Diversity/arare_max{self.seqs_min()} {self.project}/03_Diversity/Alpha'))
        os.system(f('mv {self.project}/03_Diversity/bdiv_even{self.seqs_min()} {self.project}/03_Diversity/Beta'))
        os.system(f('gzip -d {self.project}/03_Diversity/table_even{self.seqs_min()}.biom.gz'))
        os.system(f('mv {self.project}/03_Diversity/table_even{self.seqs_min()}.biom {self.project}/03_Diversity/otu_table_even.biom'))

        os.system(f('alpha_diversity.py \
                    -i {self.project}/03_Diversity/otu_table_even.biom \
                    -o {self.project}/03_Diversity/Alpha/alpha_diversity_index.txt \
                    -t {self.project}/02_OTU/otus.tre \
                    -m shannon,simpson,ace,goods_coverage,chao1,observed_species,PD_whole_tree'))

        # 生成R所需文件
        for file in os.listdir(f('{self.project}/03_Diversity/Alpha/alpha_div_collated')):
            newfile = f('{self.project}/03_Diversity/Alpha/alpha_div_collated/{file}')
            redo_data(self.project, newfile)

    def beta_diversity(self):
        # 计算贝塔多样性距离指数
        os.system(f("sed -i 's/SampleID/#SampleID/' {self.project}/sample_info.txt"))
        self.if_exists(f('{self.project}/03_Diversity/Beta/Beta_div'))
        for file in glob.glob(f('{self.project}/03_Diversity/Beta/*_dm.txt')):
            distance = '_'.join(re.split('[/_]', file)[-3:-1])
            os.system(f('make_distance_boxplots.py -d {file} \
                        -m {self.project}/sample_info.txt -f Group \
                        -o {self.project}/03_Diversity/Beta/Beta_div/{distance} \
                        --save_raw_data --suppress_all_with --suppress_all_between --suppress_individual_between'))
        os.system(f("sed -i 's/#//' {self.project}/sample_info.txt"))

    def bak(self):
        os.system(f('biom convert -i {self.project}/03_Diversity/otu_table_even.biom -o {self.project}/03_Diversity/otu_table_even.txt \
                    --table-type="OTU table" --to-tsv'))
        os.system(f('biom convert -i {self.project}/03_Diversity/otu_table_even.biom -o {self.project}/03_Diversity/otu_table_even_tax.txt \
                    --table-type="OTU table" --header-key taxonomy --to-tsv'))
        os.system(f("sed -i.bak 's/#//' {self.project}/02_OTU/otu_table.txt"))
        os.system(f("sed -i.bak '/# Const/d;s/#//' {self.project}/02_OTU/otu_table_tax.txt"))     
        os.system(f("sed -i.bak '/# Const/d;s/#//g' {self.project}/03_Diversity/otu_table_even.txt"))
        os.system(f("sed -i.bak '/# Const/d;s/#//' {self.project}/03_Diversity/otu_table_even_tax.txt"))
        #处理物种注释结果文件
        tax(self.project)


    def lefse(self):
        os.system(f("sed -i 's/SampleID/#SampleID/' {self.project}/sample_info.txt"))
        os.system(f('summarize_taxa.py -i {self.project}/03_Diversity/otu_table_even.biom \
                    -o {self.project}/03_Diversity/Beta/lefse/tmp -m {self.project}/sample_info.txt --delimiter "|"'))
        os.system(f("sed -i 's/#//' {self.project}/sample_info.txt"))
        os.system(f("sed -i 's/#*//' {self.project}/03_Diversity/Beta/lefse/tmp/*.txt"))
        ## 获取比较组
        sample_info = pd.read_csv(f('{self.project}/sample_info.txt'), sep='\t')
        groups = set(sample_info['Group'].tolist())
        groupvs_list = combinations(groups, 2)
        ## 对每个比较组分析
        for file in glob.glob(f('{self.project}/03_Diversity/Beta/lefse/tmp/*.txt')):
            if re.search('L6', file):
                L6 = pd.read_csv(file, index_col=0, sep='\t')
                for groupvs in groupvs_list:
                    tmp_groupvs = '_vs_'.join(groupvs)
                    self.if_exists(f('{self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}'))
                    tmp_L6 = L6[L6['Group'].isin(list(groupvs))]
                    tmp_file = f('{self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}/{tmp_groupvs}_L6.txt')
                    tmp_L6.to_csv(tmp_file, sep='\t')
                    out = f('{self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}/{tmp_groupvs}.txt')
                    os.system(f("awk -F'\t' '{{$2=null;$3=null;$5=null;print $0}}' {tmp_file} > {out}"))   #### 删除不必要的列
                    os.system(f("sed -i 's/\s\+/\t/g' {out}"))    ### 替换连续空格为tab
                    os.system(f('lefse-format_input.py {self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}/{tmp_groupvs}.txt {self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}/{tmp_groupvs}.in -f c -c 2 -s -1 -u 1 -o 1000000'))
                    os.system(f('run_lefse.py {self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}/{tmp_groupvs}.in {self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}/{tmp_groupvs}.res'))
                    os.system(f('lefse-plot_res.py {self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}/{tmp_groupvs}.res {self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}/{tmp_groupvs}.png --dpi 600'))
                    os.system(f('lefse-plot_cladogram.py {self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}/{tmp_groupvs}.res {self.project}/03_Diversity/Beta/lefse/{tmp_groupvs}/{tmp_groupvs}.cla.png --format png --dpi 600'))

    def main(self):
        # self.flash()
        # self.qiime_qc()
        # self.data_statis()
        # self.derep_fa()
        # self.otu()
        # self.annotation()
        # self.build_tree()
        # self.diversity()
        # self.alpha_diversity()
        self.beta_diversity()
        # self.bak()
        self.lefse()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='16s pipeline')
    parser.add_argument('-i', '--input', help='input the path of project', required=True)
    args = parser.parse_args()
    project = args.input
    project = project.rstrip('/')
    pipeline = Pipeline(project)
    pipeline.main()
    # otu.main(project)