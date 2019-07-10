import os
import sys
import re
import pandas as pd
import numpy as np
from ww import f

pd.set_option('display.width', None)
pd.set_option('display.max_columns', None)    # 设置显示所有列
pd.set_option('display.max_rows', None)      # 显示所有行

def sample_info(path):
    sample_df = pd.read_csv(f('{path}/sample_info.txt'), sep='\t')
    sample = sample_df[['SampleID', 'Group']]
    sample = sample.groupby('Group')
    group = set(sample_df['Group'].tolist())
    return sample, group

def if_exists(file):
    if not os.path.exists(file):
        os.makedirs(file)

def clean(file):
    data = pd.read_csv(file, index_col=0, sep='\t')
    # 拆分分类信息
    tax_df = data['taxonomy'].str.split(';', expand=True)
    tax_df.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    tax_df.replace('.*__$', ' Others', regex=True, inplace=True)
    tax_df = tax_df.fillna(' Others')
    # 按列合并数据（横向拼接）
    data = data.drop('taxonomy', axis=1)
    df = pd.concat([data, tax_df], axis=1)
    #对nan填充字符
    # df = df.fillna('Others')
    # print(df.isnull().any())
    return df

def classified_stat(path, df):
    # print(df.head())
    if_exists(f('{path}/02_OTU/taxa_stat'))
    index_list = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    stat_df = pd.DataFrame(columns=index_list)
    for index in index_list:
        df_group_sum = df.groupby(index).sum()
        if ' Others' in df_group_sum.index:
            df_group_sum = df_group_sum.drop(' Others')
            stat_df[index] = df_group_sum.sum()
        else:
            stat_df[index] = df_group_sum.sum()
    # print(stat_df)
    stat_df.to_csv(f('{path}/02_OTU/taxa_stat/taxa_stat.csv'))

def absolute(path, df, index):
    # 按分类类别进行分组汇总
    df_group = df.groupby(index).sum()
    others = pd.DataFrame(df_group.loc[' Others']).T
    df_group = df_group.drop(' Others')
    df_abs = pd.concat([df_group, others], axis=0)
    # print(df_group)
    if_exists(f('{path}/02_OTU/Absolute'))
    df_abs.to_csv(f('{path}/02_OTU/Absolute/{index}.csv'))

    return df_abs

def relative(path, df, index):
    # 计算分类后的每个样本含有类别的总数
    sample_total = pd.DataFrame(df.sum()).T
    sample_total.rename(index={0: 'Total'}, inplace=True)
    df_tmp = pd.concat([df, sample_total])
    # 计算相对值
    df_rela = df_tmp.apply(lambda x: x / x['Total'], axis=0)
    df_rela = df_rela.drop('Total', axis=0)
    if_exists(f('{path}/02_OTU/Relative'))
    df_rela.to_csv(f('{path}/02_OTU/Relative/{index}.csv'))

    return df_rela

    # 聚类热图数据
def cluster(path, df, index):
    if_exists(f('{path}/02_OTU/taxa_heatmap/cluster'))
    if index != 'Species':
        tax_cluster = df.drop(' Others')
        tax_cluster.to_csv(f('{path}/02_OTU/taxa_heatmap/cluster/{index}.csv'))

def top10(path, df, index):
    #对每一个类进行汇总（求行和）并排序
    tax_sort = df.sum(axis=1).sort_values(ascending=False)
    # top10的类名称及其他
    tax_list = tax_sort.index.tolist()
    if ' Others' not in tax_list[0:10]:
        top10 = tax_list[0:10]
        other = tax_list[10:]
    else:
        tax_list.remove(' Others')
        top10 = tax_list[0:10]
        other = tax_list[10:].append(' Others')
    # 按index取出top10和其他
    top10_df = df.loc[top10].T
    if other:
        top10_df[' Others'] = df.loc[other].sum()
    # print(top10_df)
    # 添加total
    top10_df['Total'] = top10_df.sum(axis=1)
    # 求相对值
    tax10_df = top10_df.apply(lambda x: x / x['Total'], axis=1)
    tax10_df = tax10_df.drop('Total', axis=1)
    # print(tax10_df)
    if_exists(f('{path}/02_OTU/Top10'))
    tax10_df.to_csv(f('{path}/02_OTU/Top10/{index}_top10.csv'))

    return tax_list, tax10_df

def group(path, df1, df2, index, tax_list):
    # 按组计算
    tax_group = pd.DataFrame(index=tax_list, columns=list(sample_info(path)[1]))
    tax10_group = pd.DataFrame(columns=top10)
    for group in sample_info(path)[0]:
        sample_name = group[1]['SampleID'].tolist()
        # 对组内求均值
        tax_tmp = pd.DataFrame(df1[sample_name].mean(axis=1))
        tax_tmp.rename(columns={0: group[0]}, inplace=True)
        tax_group[group[0]] = tax_tmp[group[0]]
        tax10_tmp = pd.DataFrame(df2.loc[sample_name].mean()).T
        tax10_tmp.rename(index={0: group[0]}, inplace=True)
        tax10_group = pd.concat([tax10_group, tax10_tmp], sort=False)
    tax_group.sort_index(inplace=True)
    # print(tax_group)
    if_exists(f('{path}/02_OTU/Abundance/Relative_group'))
    tax_group.to_csv(f('{path}/02_OTU/Abundance/Relative_group/{index}_group.csv'))
    if_exists(f('{path}/02_OTU/taxa_heatmap/cluster_group'))
    if 'Others' in tax_group.index:
        tax_group_cluster = tax_group.drop('Others')
        tax_group_cluster.to_csv(f('{path}/02_OTU/taxa_heatmap/cluster_group/{index}_group.csv'))
    else:
        tax_group.to_csv(f('{path}/02_OTU/taxa_heatmap/cluster_group/{index}_group.csv'))
    if_exists(f('{path}/02_OTU/Top10_group'))
    tax10_group.to_csv(f('{path}/02_OTU/Top10_group/{index}_group_top10.csv'))


def write(path, df, index, num):
    # 按分类类别进行分组汇总
    df_group = df.groupby(index).sum()
    others = pd.DataFrame(df_group.loc[' Others']).T
    df_group = df_group.drop(' Others')
    df_group = pd.concat([df_group, others], axis=0)
    # print(df_group)
    if_exists(f('{path}/02_OTU/Abundance/Absolute'))
    df_group.to_csv(f('{path}/02_OTU/Abundance/Absolute/{index}.csv'))
    # 添加分类详细信息
    # tax = df['Kingdom'].str.cat([df['Phylum'], df['Class'],
    #                             df['Order'], df['Family'],
    #                             df['Genus'], df['Species']], sep=';')
    # tax = tax.drop_duplicates()
    # tmp_index = [i for i in df_group.index]
    # tmp_index.remove(' Others')
    # tax_detail_list = {';'.join(i.split(';')[:num]) for i in tax}
    # tax_detail = [j for i in tmp_index for j in tax_detail_list if i == j.split(';')[-1]]
    # tax_detail = sorted(set(tax_detail), key=tax_detail.index)
    # tax_detail.append(' Others')
    # print(len(tmp_index))
    # print(len(tax_detail))
    # if len(df_group.index) == len(tax_detail):
    #     df_group['tax_detail'] = tax_detail
    # if_exists(f('{path}/02_OTU/Absolute'))
    # df_group.to_csv(f('{path}/02_OTU/Absolute/{index}.csv'))
    # if 'tax_detail' in df_group.columns:
    #     df_group = df_group.drop('tax_detail', axis=1)

    # print(df_group)
    # 计算分类后的每个样本含有类别的总数
    sample_total = pd.DataFrame(df_group.sum()).T
    sample_total.rename(index={0: 'Total'}, inplace=True)
    tax_df = pd.concat([df_group, sample_total])
    # 计算相对值
    tax_df = tax_df.apply(lambda x: x / x['Total'], axis=0)
    tax_df = tax_df.drop('Total', axis=0)
    if_exists(f('{path}/02_OTU/Abundance/Relative'))

    tax_df.to_csv(f('{path}/02_OTU/Abundance/Relative/{index}.csv'))
    # 聚类热图数据
    if_exists(f('{path}/02_OTU/taxa_heatmap/cluster'))
    if index != 'Species':
        tax_cluster = tax_df.drop(' Others')
        # if len(tax_cluster.index) == len(tax_detail):
        #     tax_cluster['tax_detail'] = tax_detail
        tax_cluster.to_csv(f('{path}/02_OTU/taxa_heatmap/cluster/{index}.csv'))

    #对每一个类进行汇总（求行和）并排序
    tax_sort = df_group.sum(axis=1).sort_values(ascending=False)
    # top10的类名称及其他
    tax_list = tax_sort.index.tolist()
    if ' Others' not in tax_list[0:10]:
        top10 = tax_list[0:10]
        other = tax_list[10:]
    else:
        tax_list.remove(' Others')
        top10 = tax_list[0:10]
        other = tax_list[10:].append(' Others')
    # 按index取出top10和其他
    top10_df = df_group.loc[top10].T
    if other:
        top10_df[' Others'] = df_group.loc[other].sum()
    # print(top10_df)
    # 添加total
    top10_df['Total'] = top10_df.sum(axis=1)
    # 求相对值
    tax10_df = top10_df.apply(lambda x: x / x['Total'], axis=1)
    tax10_df = tax10_df.drop('Total', axis=1)
    # print(tax10_df)
    # os.makedirs(f'C:/Users/jbwang/Desktop/OTU/Top10')
    if_exists(f('{path}/02_OTU/Top10/{index}'))
    tax10_df.to_csv(f('{path}/02_OTU/Top10/{index}/{index}_top10.csv'))

    # 按组计算
    tax_group = pd.DataFrame(index=tax_list, columns=list(sample_info(path)[1]))
    tax10_group = pd.DataFrame(columns=top10)
    for group in sample_info(path)[0]:
        sample_name = group[1]['SampleID'].tolist()
        # 对组内求均值
        tax_tmp = pd.DataFrame(tax_df[sample_name].mean(axis=1))
        tax_tmp.rename(columns={0: group[0]}, inplace=True)
        tax_group[group[0]] = tax_tmp[group[0]]
        tax10_tmp = pd.DataFrame(tax10_df.loc[sample_name].mean()).T
        tax10_tmp.rename(index={0: group[0]}, inplace=True)
        tax10_group = pd.concat([tax10_group, tax10_tmp], sort=False)
    tax_group.sort_index(inplace=True)
    # print(tax_group)
    if_exists(f('{path}/02_OTU/Abundance/Relative_group'))
    tax_group.to_csv(f('{path}/02_OTU/Abundance/Relative_group/{index}_group.csv'))
    if_exists(f('{path}/02_OTU/taxa_heatmap/cluster_group'))
    if 'Others' in tax_group.index:
        tax_group = tax_group.drop('Others')
        # if len(tax_group.index) == len(tax_detail):
        #     tax_group['tax_detail'] = tax_detail
        tax_group.to_csv(f('{path}/02_OTU/taxa_heatmap/cluster_group/{index}_group.csv'))
    else:
        # if len(tax_group.index) == len(tax_detail):
        #     tax_group['tax_detail'] = tax_detail
        tax_group.to_csv(f('{path}/02_OTU/taxa_heatmap/cluster_group/{index}_group.csv'))
    if_exists(f('{path}/02_OTU/Top10_group/{index}'))
    tax10_group.to_csv(f('{path}/02_OTU/Top10_group/{index}/{index}_group_top10.csv'))

def main(path):    
    file = path + '/02_OTU/otu_table_tax.txt'
    df = clean(file)
    classified_stat(path, df)
    index_list = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    for num, index in enumerate(index_list):
        num += 2
        write(path, df, index, num)
        # df_absolute = absolute(path, df, index)
        # df_relative = relative(path, df_absolute, index)
        # cluster(path, df_relative, index)
        # tax_list = top10(path, df_absolute, index)[0]
        # tax10_df = top10(path, df_absolute, index)[1]
        # group(path, df_relative, tax10_df, index, tax_list)


if __name__=='__main__':
    path = sys.argv[1].rstrip('/')
    main(path)

