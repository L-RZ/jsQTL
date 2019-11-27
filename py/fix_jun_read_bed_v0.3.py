# This script is used to merge the junction.bed file into one.
# And get the junction position within two exons boundary
# And merge the junction read from each sample into a big table

# v0.2 remove the strand + or -
# the strand of junction will detected in exon-skippping step
# output each junction with two line of strand + and -

# v0.3 change the output method to write each line instead of a big list for reduce the memory

import os, sys


def get_intron_border(bed_line_l):
    jun_start = int(bed_line_l[1]) + int(bed_line_l[10].split(',')[0])    #start pos in exon1
    jun_end = int(bed_line_l[1]) + int(bed_line_l[11].split(',')[1]) + 1  #start pos in exon2
    return jun_start, jun_end


def make_junct_dict(bed_in_addr, junct_dict, sample_id):
    f_bed = open(bed_in_addr)
    for each_line in f_bed:
        each_line_l = each_line.split('\t')
        j_chr = each_line_l[0]
        j_start, j_end = get_intron_border(each_line_l)
        # j_strand = each_line_l[5]
        j_strand = '.'
        j_reads = int(each_line_l[4])   # 2018-9-5 fixed a bug
        j_start_bin = int(j_start / 1e6)
        j_s_e = '_'.join([str(j_start), str(j_end)])
        # if j_s_e == '155032819_155033239': # for test
            # pass
        if j_chr in junct_dict[j_strand]:
            if j_start_bin in junct_dict[j_strand][j_chr]:
                if j_s_e in junct_dict[j_strand][j_chr][j_start_bin]:
                    if sample_id in junct_dict[j_strand][j_chr][j_start_bin][j_s_e]:
                        junct_dict[j_strand][j_chr][j_start_bin][j_s_e][sample_id] += j_reads  # 2018-9-5 fixed a bug
                    else:
                        junct_dict[j_strand][j_chr][j_start_bin][j_s_e][sample_id] = j_reads
                else:
                    junct_dict[j_strand][j_chr][j_start_bin][j_s_e] = {sample_id: j_reads}
            else:
                # junct_dict[j_strand][j_chr][j_start_bin] = {}
                # junct_dict[j_strand][j_chr][j_start_bin][j_s_e] = {sample_id: j_reads}
                junct_dict[j_strand][j_chr][j_start_bin] = {j_s_e: {sample_id: j_reads}}
        else:
            junct_dict[j_strand][j_chr] = {j_start_bin: {j_s_e: {sample_id: j_reads}}}
    return junct_dict


def scan_each_bed(bed_in_addr_list, out_addr):
    n_sample = len(bed_in_addr_list)
    bed_id = 0   # bed file index
    # all_junct_dict= {'-': {},'+': {}}
    all_junct_dict = {'.': {}}
    sample_names = []
    for each_bed_addr in bed_in_addr_list:
        all_junct_dict = make_junct_dict(each_bed_addr, all_junct_dict, bed_id)
        bed_id += 1

        path, name = os.path.split(each_bed_addr)
        name, ext = os.path.splitext(name)
        sample_names.append(name)
    # output
    sample_names_str = '\t'.join(sample_names)
    header = 'Chr\tstart\tend\tstrand\t%s'%sample_names_str
    f_out = open(out_addr, 'w')
    f_out.write(header + '\n')
    # result = [header]

    for each_strand in all_junct_dict:
        strand_dict = all_junct_dict[each_strand]
        for each_chr in strand_dict:
            chr_dict = strand_dict[each_chr]
            for each_bin in chr_dict:
                bin_dict = chr_dict[each_bin]
                for each_junct in bin_dict:
                    junct_dict = bin_dict[each_junct]
                    read_list = [0] * n_sample
                    for each_sample_index in junct_dict:
                        read_list[each_sample_index] = junct_dict[each_sample_index]
                    read_str = '\t'.join(map(str, read_list))
                    each_junct_pos = each_junct.split('_')
                    # out_line = "{0}\t{1}\t{2}\t{3}\t{4}".format(each_chr, each_junct_pos[0],each_junct_pos[1],each_strand,read_str)
                    out_line_pos = "{0}\t{1}\t{2}\t{3}\t{4}".format(each_chr, each_junct_pos[0], each_junct_pos[1], '+', read_str)
                    out_line_neg = "{0}\t{1}\t{2}\t{3}\t{4}".format(each_chr, each_junct_pos[0], each_junct_pos[1], '-', read_str)

                    f_out.write(out_line_pos + '\n')
                    f_out.write(out_line_neg + '\n')

#test()
def usage():
    print '''
    e.g.
    %s output.txt input1.bed input2.bed ...'''%sys.argv[0]


def main():
    if sys.argv < 3:
        usage()
        sys.exit()
    file_in_addr_list = sys.argv[2:] # *.bed 
    file_out_addr = sys.argv[1]      
    scan_each_bed(file_in_addr_list, file_out_addr)

if __name__ == '__main__':
    main()





