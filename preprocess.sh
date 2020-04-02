#!/usr/bin/env bash

set -e;

project=$1
project=${project%/}   # 去掉斜杠
ori_path=$(pwd)
[ $project = "." ] && project=$(pwd)
cd $project
source activate qiime1
## 拆分数据
forward_primer=$(sed -n '2p' sample_info.txt | cut -f3)
reverse_primer=$(sed -n '2p' sample_info.txt | cut -f4)
reverse_primer2=$(echo $reverse_primer | tr 'ATCGRYMKBVDH' 'TAGCYRKMVBHD' | rev) # 取反向互补
# echo $reverse_primer2

for i in $(ls ${project}/00_RawData/*/* | grep '1.fq.gz')
do 
	# echo $i
	j=${i/1.fq.gz/2.fq.gz}
	# echo $j
	fq=${i##*/}
	tmp=$(echo $fq | cut -d "-" -f1)
	sample_group=$(echo $fq | cut -d "_" -f1)
	sample_group=${sample_group/${tmp}-/}
	echo $sample_group
	mkdir -p 00_RawData/Extend/${sample_group}_split
	/home/jbwang/soft/FLASH-1.2.11/flash $i $j -d 00_RawData/Extend -o $sample_group -x 0.1 -M 150 -z && \
	/home/jbwang/soft/fastq-multx/fastq-multx -m 1 -B ${sample_group}-barcode.txt -b 00_RawData/Extend/${sample_group}.extendedFrags.fastq.gz -o 00_RawData/Extend/${sample_group}_split/%.extendedFrags.fastq.gz
	
	# # 去引物

	rm -rf 00_RawData/Extend/${sample_group}_split/unmatched*
	mkdir -p 00_RawData/Extend/${sample_group}_noprimer/tmp

	# # echo 'GACTACHVGGGTATCTAATCC' | tr 'ATCGRYMKBVDH' 'TAGCYRKMVBHD' | rev  
	for s in $(ls 00_RawData/Extend/${sample_group}_split)
	do
		# echo $s
		cutadapt -g $forward_primer -e 0.15 00_RawData/Extend/${sample_group}_split/$s -o 00_RawData/Extend/${sample_group}_noprimer/tmp/$s && \
		cutadapt -a $reverse_primer2 -e 0.15 00_RawData/Extend/${sample_group}_noprimer/tmp/$s -o 00_RawData/Extend/${sample_group}_noprimer/$s
	done

	# # 质控
	mkdir -p 00_RawData/TrimQC/tmp
	for q in $(ls 00_RawData/Extend/${sample_group}_noprimer | grep 'gz')
	do
		sample=$(echo $q | cut -d "." -f1)
		# echo $sample
		split_libraries_fastq.py  -i 00_RawData/Extend/${sample_group}_noprimer/${q} --sample_ids ${sample} -o 00_RawData/TrimQC/${sample} \
                        -q 19 --max_bad_run_length 3 --min_per_read_length_fraction 0.75 --max_barcode_errors 0 \
                        --store_demultiplexed_fastq --barcode_type not-barcoded --phred_offset 33 && \
        mv 00_RawData/TrimQC/${sample}/seqs.fna 00_RawData/TrimQC/${sample}/${sample}.fa
        mv 00_RawData/TrimQC/${sample}/seqs.fastq 00_RawData/TrimQC/${sample}/${sample}.fq

        # 对每个样品去重
        /home/jbwang/soft/Usearch/Usearch11 -fastx_uniques 00_RawData/TrimQC/${sample}/${sample}.fa -fastaout 00_RawData/TrimQC/${sample}/${sample}_norep.fa && \
        # Usearch检测嵌合体
        /home/jbwang/soft/Usearch/Usearch11 -uchime2_ref 00_RawData/TrimQC/${sample}/${sample}_norep.fa -db /home/jbwang/refedata/gold/gold.fa \
                        -chimeras 00_RawData/TrimQC/tmp/chimeras.fa -notmatched 00_RawData/TrimQC/tmp/notgold.fa -uchimeout 00_RawData/TrimQC/tmp/gold.uchime \
                        -strand plus -mode balanced -threads 16 && \
        # 获得嵌合体序列ID
        grep '>' 00_RawData/TrimQC/tmp/chimeras.fa| sed 's/>//g' > 00_RawData/TrimQC/tmp/chimeras.id
        mkdir -p 01_CleanData/${sample}
        filter_fasta.py -f 00_RawData/TrimQC/${sample}/${sample}.fa -o 01_CleanData/${sample}/${sample}.fa -s 00_RawData/TrimQC/tmp/chimeras.id -n && \
        sed -i 's/ .*//g' 01_CleanData/${sample}/${sample}.fa
        # 获得所有序列的ID
        grep '>' 00_RawData/TrimQC/${sample}/${sample}.fa | sed 's/>//g' > 00_RawData/TrimQC/tmp/ids
        # 去除嵌合体的ID
        grep -F -v -f 00_RawData/TrimQC/tmp/chimeras.id 00_RawData/TrimQC/tmp/ids | sed 's/>//g' > 00_RawData/TrimQC/tmp/id
        # 筛选fastq
        /home/jbwang/soft/seqtk/seqtk subseq 00_RawData/TrimQC/${sample}/${sample}.fq 00_RawData/TrimQC/tmp/id > 01_CleanData/${sample}/${sample}.fq && \
        # rm 00_RawData/TrimQC/tmp/chimeras.id 00_RawData/TrimQC/tmp/chimeras.fa 00_RawData/TrimQC/tmp/notgold.fa 00_RawData/TrimQC/tmp/gold.uchime 00_RawData/TrimQC/tmp/ids 00_RawData/TrimQC/tmp/id
        rm 00_RawData/TrimQC/tmp/*
        rm -rf 01_CleanData/${sample}/sed*
	done
done

python3 ~/code/QcStatic.py .
cd $ori_path

