import os
import sys
import glob
import q30
from q30 import read_count
from q30 import stat

def get_fq_file(project):
    rawPath = os.path.join(project, '00_RawData', 'Extend')
    cleanPath = os.path.join(project, '01_CleanData')
    rawFile = [file for file in glob.glob(rawPath + '/*_split/*.gz')]
    mydict = {i.split('/')[-1]:'/'.join(i.split('/')[:-1]) for i in rawFile}
    rawFile_list = sorted([i.split('/')[-1] for i in rawFile])
    rawFile_list = [os.path.join(mydict.get(fq), fq) for fq in rawFile_list]
    # rawFile_list = sorted([project + '/00_RawData/' + file for file in os.listdir(rawPath) if re.search(r'_?R?1.f(ast)?q.?(gz)?', file)])
    cleanFile_list = sorted([file for file in glob.glob(cleanPath + '/*/*.fq')])

    return rawFile_list, cleanFile_list

def data_statis(project):
    print('QC Static...')
    static_file = project + '/QcStatic.csv'
    with open(static_file, 'w') as f:
        f.write('SampleID,Raw PE,Clean PE,Base(nt),AvgLen(nt),Q20(%),Q30(%),GC(%),Effective(%)\n')
        for rawFile, cleanFile in zip(get_fq_file(project)[0], get_fq_file(project)[1]):
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
    print('Static Done')

if __name__ == '__main__':
	project = sys.argv[1].rstrip('/')
	data_statis(project)

