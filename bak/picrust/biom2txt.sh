#!/bin/bash

biom=$1

## 替换文件名后缀biom为txt
out=(${biom//biom/txt})
#echo $out
biom convert -i $biom -o $out --table-type="OTU table" --to-tsv
