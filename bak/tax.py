import sys
from ww import f

def tax(path):
    with open(f('{path}/02_OTU/otus_tax_assignments.txt')) as f1:
        with open(f('{path}/02_OTU/tax.txt'), 'w') as f2:
            # f2.write(f1.readline())
            for line in f1:
                otu_id = line.split('\t')[0]
                tax = line.strip().split('\t')[1].split(';')
                if len(tax) == 7:
                    f2.write(otu_id + '\t' + '\t'.join(tax) + '\n')
                else:
                    while len(tax) < 7:
                        tax.append('NA')
                    newline = otu_id + '\t' + '\t'.join(tax)
                    f2.write(newline + '\n')
                if len(tax) > 7:
                    print(line)

path = sys.argv[1].rstrip('/')
tax(path)
