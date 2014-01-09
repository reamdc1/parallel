from multiprocessing import Pool
import numpy
import time
import sys
import os
import Bio
from Bio import SeqIO,SeqFeature
from Bio.SeqRecord import SeqRecord


# So the purpose of this bit of code is to determine possible targets of HGT for the first paper.  This is done by generating a summary 
# by making CSV files wich have the local operon structure. The resulst are organized in the order of phylogenetic distance. We will suspect
# that there is HGT if an operon is found in one organism whose branch does not contain the operon. 

#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result

def return_filtered_information(fname):
    print fname
    result = {}
    for line in [i.strip() for i in open(fname).readlines()]:
        tmp = line.split('\t')
        if tmp[0] in ['@@', '$$', '##']:
            pass
        elif tmp[0] == '++':
            result.update({nc: "%s,%s,%s" % (nc, org_name, tmp[1]) })
        elif len(tmp) > 1:
            #print tmp
            nc = tmp[0]
            org_name = tmp[1]
            
    return result
            

def main():
    #print "hi there"
    folder = './optimized_results_proteobacteria'
    result_folder = './operon_structure_result/'
    
    order_list = [i.strip() for i in open('./phylo_order.txt').readlines()]
    #print order_list
    flist = returnRecursiveDirFiles(folder)
    for fname in flist:
        print fname
        operon = fname.split('/')[len(fname.split('/'))-1].split('.')[0]
        print operon
        res = return_filtered_information(fname)
        result =  [operon]
        for nc in order_list:
            if nc in res.keys():
                result.append(res[nc])
        print result
        handle = open(result_folder + operon + '.csv', 'w')
        handle.write('\n'.join(result))
        handle.close()
            
        














if __name__ == '__main__':
    main()
