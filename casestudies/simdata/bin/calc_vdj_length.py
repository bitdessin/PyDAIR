import sys
import os
import re
import numpy as np

v_file = sys.argv[1]
d_file = sys.argv[2]
j_file = sys.argv[3]


def calc_len(f):
    gene_len = []
    with open(f, 'r') as fh:
        for buf in fh:
            if buf[0:1] == '>':
                m = re.search("([0-9]+) nt", buf)
                if m:
                    gene_len.append(int(m.group(1)))
                else:
                    raise ValueError('missed in sequence length?')
    gene_len = np.array(gene_len)
    return gene_len

v_len = calc_len(v_file)
d_len = calc_len(d_file)
j_len = calc_len(j_file)

print "Average length of V gene:" + str(np.mean(v_len)) + " std: " + str(np.std(v_len)) 
print "Average length of D gene:" + str(np.mean(d_len)) + " std: " + str(np.std(d_len))
print "Average length of J gene:" + str(np.mean(j_len)) + " std: " + str(np.std(j_len))


