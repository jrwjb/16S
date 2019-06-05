import os
from ww import f

def otu2biom(path, otu_table):
    # otu_table.txt 转换为 biom 格式
    biom = otu_table.replace('txt', 'biom')
    biom_tax = otu_table.split('.')[0] + '_tax.biom'
    otu_table_tax = biom_tax.replace('biom', 'txt')
    sum = biom_tax.replace('biom', 'sum')
    os.system(f('biom convert -i {path}/{otu_table} \
                        -o {path}/{biom} \
                        --table-type="OTU table" \
                        --to-json'))

    # # 添加物种信息至OTU表最后一列，命名为taxonomy
    os.system(f('biom add-metadata -i {path}/{biom} \
                        --observation-metadata-fp {path}/otus_tax_assignments.txt \
                        -o {path}/{biom_tax} \
                        --sc-separated taxonomy \
                        --observation-header OtuId,taxonomy'))

    os.system(f('biom convert -i {path}/{biom_tax} \
                        -o {path}/{otu_table_tax} \
                        --table-type="OTU table" \
                        --header-key taxonomy \
                        --to-tsv'))

    # # 查看OTU表的基本信息：样品，OTU数量统计
    os.system(f('biom summarize-table \
                -i {path}/{biom_tax} \
                -o {path}/{sum}'))