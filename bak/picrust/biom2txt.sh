#!/bin/bash

biom=$1

## �滻�ļ�����׺biomΪtxt
out=(${biom//biom/txt})
#echo $out
biom convert -i $biom -o $out --table-type="OTU table" --to-tsv
