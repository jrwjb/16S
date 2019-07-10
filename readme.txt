16s pipeline

1、激活环境
source activate qiime1

2、执行脚本
python3 main.py -i project -db g

参数说明：-i 输入为项目路径 -db 为参考数据库，g为greengenes, s为silva，默认为greengene。如果是默认数据库，可不输入-db参数。
另：可执行 python3 main.py -h 查看帮助。

3、绘图
需先退出环境
1）conda deactivate
2）python3 plot.py project   ### project参数为项目路径

