# to find the min pvalue SNP in exon-skippping and output the that line
import sys
print ' '.join(sys.argv)
in_addr = sys.argv[1] # *geno_gp_p_sorted_uniq.txt
out_addr = sys.argv[2]
f_in  = open(in_addr)
f_out = open(out_addr, 'w')
gene_id_tmp = ''
header = f_in.next()
f_out.write(header)
for line in f_in:
    l_line = line.strip().split()
    gene_id = l_line[0]
    # p_val = float(l_line[5])  # for whole blood np col
    p_val = float(l_line[4])    # for np only test result
    if gene_id_tmp == '':
        min_p = p_val
        min_line = line
        gene_id_tmp = gene_id
    elif gene_id != gene_id_tmp:
        f_out.write(min_line)
        gene_id_tmp = gene_id
        min_line = line
        min_p = p_val
    else:
        if min_p > p_val:
            min_p = p_val
            min_line = line
else:
    f_out.write(min_line)



