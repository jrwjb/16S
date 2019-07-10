import re
import pandas as pd
from ww import f

def redo_data(path, file):
    def sample_info(path):
        df = pd.read_csv(path + '/sample_info.txt', sep='\t')
        sample_info_dict = {i:j for i,j in zip(df['SampleID'], df['Group'])}
        return sample_info_dict

    index = re.split(r'[/.]', file)[-2]
    newfile = '/'.join(re.split(r'[/.]', file)[:-1]) + '_T.txt'
    # print(newfile)
    with open(file) as f1:             
        with open(newfile, 'w') as f2:       
            f2.write(f('sampleID\tGroup\tSeqs\t{index}\n'))
            colnames = f1.readline()
            samplenames = colnames.strip('').split()[4:]

            for line in f1:
                line = line.strip().split()
                Seqs = line[1]
                index_list = line[3:]
                for i, j in enumerate(samplenames):
                    # group = re.search('[a-zA-Z]*', j).group()
                    group = sample_info(path)[j]
                    index_data = index_list[i]
                    f2.write(f('{j}\t{group}\t{Seqs}\t{index_data}\n'))

if __name__ == '__main__':
    file = '/home/apt/Desktop/chao1.txt'
    redo_data(file)