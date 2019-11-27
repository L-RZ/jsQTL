# extract sig GeneID based on the sig count table

import sys, os
in_sig_uniq_addr = sys.argv[1]   # sig count table
in_p_addr = sys.argv[2]
name, ext = os.path.splitext(in_p_addr)

out_addr = name + '_sig.txt'

f_in_sig = open(in_sig_uniq_addr)
sig_list = [line.split()[0] for line in f_in_sig]
sig_list_set = set(sig_list)

f_in_p = open(in_p_addr)
f_out = open(out_addr, 'w')
for line in f_in_p:
    if line.split()[0] in sig_list_set:
        f_out.write(line)

