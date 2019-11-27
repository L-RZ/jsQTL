## This script is for splitting the merge_junction_reads_counts by Chr to increase the speed
# split_junct_by_chr.py junct_extract_merge.sorted.txt
# will output junct_extract_merge.sorted_1.txt ....
import os,sys

def split_file(in_addr):
    f_in = open(in_addr)
    dat = f_in.readlines()
    # chr_id_list = []
    header = dat[0]
    # result = [header]
    chr_id_old = ''
    for each in dat[1:]:
        each_l = each.split('\t')
        chr_id = each_l[0]
        if chr_id_old != '':
            if chr_id == chr_id_old:
                result.append(each)
            else:
                name, ext = os.path.splitext(in_addr)
                out_addr = '%s_%s.txt'%(name,chr_id_old)
                f_out = open(out_addr,'w')
                f_out.write(''.join(result))
                f_out.close()
                result = [header,each]
                chr_id_old = chr_id
        else:
            result = [header, each]
            chr_id_old = chr_id
    else:
        name, ext = os.path.splitext(in_addr)
        out_addr = '%s_%s.txt' % (name, chr_id_old)
        f_out = open(out_addr, 'w')
        f_out.write(''.join(result))
        f_out.close()

def main():
    print ' '.join(sys.argv)
    in_addr = sys.argv[1]
    split_file(in_addr)
    print 'Done'



main()
