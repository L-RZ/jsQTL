# used for filter exon-skipping event based on Num of leave-one-out pval < 0.05
# That means counting of num of sample with leave-one-out pval < 0.05 on each exon-skipping event
# filter : N >= 3


import sys
import numpy as np


in_addr = sys.argv[1] # *_filterInOut_para.txt 
out_addr = sys.argv[2]

f_in = open(in_addr)
header = f_in.next()
f_out = open(out_addr, 'w')
f_out.write(header.split()[0] + '\t' + 'N_0.05\n')
for i, line in enumerate(f_in):
    if i % 2 == 0:
        l_line = line.strip().split()
        l_line_p = l_line[3:]
        n = np.sum(np.array(l_line_p, dtype=float) < 0.05)
        f_out.write(l_line[0] + '\t' + str(n) + '\n')

