# v0.1
# remove exon-skipping event by Number of outliers which are from leave-one-out p-val < 0.05
#

import sys

in_addr = sys.argv[1]  # whole_blood_junct_exonskipping_-_1_sumLJ_filterInOut.txt
in_N_addr = sys.argv[2] # whole_blood_junct_exonskipping_-_1_sumLJ_filterInOut_para_n0.05.txt
out_addr = sys.argv[3]

N_set = 5 # N_outlier > 5

f_in_N = open(in_N_addr)
header = f_in_N.next()
N_outlier_dict = {}
for line in f_in_N:
    line_l = line.strip().split('\t')
    id_exon = line_l[0]
    N_outlier = line_l[1]
    N_outlier_dict[id_exon] = int(N_outlier)
f_in_N.close()

f_in = open(in_addr)
f_out = open(out_addr, 'w')

header = f_in.next()
f_out.write(header)

for line in f_in:
    line_l = line.strip().split('\t')
    id_exon = line_l[0]
    if id_exon in N_outlier_dict:
        N_outlier = N_outlier_dict[id_exon]
        if N_outlier > N_set:
            f_out.write(line)

f_in.close()
f_out.close()



