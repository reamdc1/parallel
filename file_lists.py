#!/usr/bin/python

# Copyright(C) 2013 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment

#from multiprocessing import Pool
#import time
import os
import simplejson as json
import argparse
import gc

from multiprocessing import Pool
import time
import sys
import os
import re
import Bio
from Bio import SeqIO,SeqFeature
from Bio.SeqRecord import SeqRecord
import itertools
import math
from Bio.SeqUtils import GC




NUM_PROCESSORS = os.sysconf("SC_NPROCESSORS_CONF")


#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# This function will take in the folder where the genomes reside.  By deafult i have this under Desktop/all_genbank/
# but this could be different. This function will make a master file list, genome only file list, master organism file list,
# and a filtered set. Filtering will be at the NC number level right now. I am not allowing renaming of the file names, because
# there is really no reason for this nonsense.
def return_genbank_paths_of_interest(in_folder = '/home/dave/Desktop/all_genbank', out_folder = './genbank_pathway_lists/', filter_fname = ''):
    print "in_folder", in_folder, "out_folder", out_folder, "filter_fname", filter_fname
    
    
    flist = returnRecursiveDirFiles(in_folder)
    
    # This first block is just the raw dump of a recursive folder walk.
    handle = open(out_folder+"all_genbank_paths.txt", 'w')
    handle.write('\n'.join(flist))
    handle.close()
    
    # This block of code is to provide a list of chromosomes and plasmids that are contained within an organisms' folder
    # The formatting will be organism name (as read from the folder name, and generally/always? corresponds to the common
    # name of the organism, then a tab delineated list of file paths.
    f_dict = {}
    for fname in flist:
        folder = fname.split('/')[-2]
        if folder not in f_dict.keys():
            f_dict.update({folder:[fname]})
        else:
            f_dict[folder].append(fname)
    #print len(f_dict)
    #print f_dict.keys()
    
    handle = open(out_folder+"grouped_by_organism_all_genbank_paths.txt", 'w')
    for org in sorted(f_dict.keys()):
        lst = [org] + sorted(f_dict[org])
        handle.write('\t'.join(lst)+'\n')
    handle.close()
    
    
    handle = open(out_folder+"filtered_list_genbank_paths.txt", 'w')
    if filter_fname != '':
        filter_list = [i.strip() for i in open(filter_fname).readlines()]
        res = []
        for fname in flist:
            if fname.split('/')[-1].split('.')[0] in filter_list:
                res.append(fname)
    handle.write('\n'.join(res))
    handle.close()
    
    
def convert_genbank(genbank_path):
    #genbank_path, hgt_dict, db_directory = genbank_tuple
    #print 'hooty hoo', genbank_path#, hgt_dict, db_directory
    #record_list = []
    
    handle = open(genbank_path)
    
    #seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
    #seq_record = SeqIO.parse(handle, "genbank").next()
    #seq_record = SeqIO.parse(handle, "genbank").next()
    
    #seq_record = SeqIO.parse(genbank_path, "genbank").next()
    
    seq_record = SeqIO.parse(handle, "genbank").next()
    #del SeqIO.parse
    #print "length", len(seq_record)
    
    handle.close()
    del handle
    #gc.collect()
    nc_number = seq_record.name
    
    # genbank information about the sequence found in the file. 
    accession = seq_record.id.split('.')[0]
    comment = seq_record.annotations['comment']
    data_file_division = seq_record.annotations['data_file_division'] # this is not very useful, since it is historical and does not reflect current taxonomy
    date = seq_record.annotations['date'] # this should be the date that the file was updated last. uesful if i wish to autoupdate for more current versions of a genome
    gi = seq_record.annotations['gi'] # i do not plan on using this, but it may be handly to have, so i will include anyway
    key_words = seq_record.annotations['keywords']
    organism = seq_record.annotations['organism'].replace(' ', '_')
    sequence_version = seq_record.annotations['sequence_version'] # this is an integer
    taxonomy = seq_record.annotations['taxonomy'] # this will be a list
    source_folder = genbank_path.split('/')[len(genbank_path.split('/'))-2]
    
    # this block of code determines the type of sequence [chromosome, plasmid] and if there are more than one chromosome in the organism.
    # In the case of mutiple chromosomes, the code will assign a name to each, based on the information provided in the genbank file.
    description = seq_record.description
    del seq_record
    #print description
    if re.search("complete genome", description) is not None: # This is a complete genome
        #print 'complete_genome'
        seq_type = 'complete_genome'
        seq_name = 'complete_genome'
    elif re.search("hromosome", description) is not None: # we have a segmented chromosome (more than one chromosome)
        seq_type = 'chromosome'
        seq_name = 'chromosome' + description.split('hromosome')[1].split(',')[0].strip()
        #print 'seq_name', seq_name
    elif re.search("lasmid", description) is not None: # we have a plasmid
        seq_type = 'plasmid'
        seq_name = description.split('lasmid')[1].split(',')[0].strip()
    else:
        print 'error', accession
        print description
        seq_type = 'complete_genome'
        seq_name = 'complete_genome'


    #handle.close()
    #gc.collect()
    #del seq_record
    gc.collect()
    return accession, {'seq_type': seq_type, 'seq_name': seq_name, 'taxonomy': taxonomy, 'organism': organism, 'source_folder': source_folder, 'file_path': genbank_path}


def parallel_convert_genbank(genbank_path_list, outfile = './genbank_pathway_lists/nc_information_dict.json'):

    pool = Pool(processes = NUM_PROCESSORS)
    #pool = Pool(processes = 1)
    result = dict(pool.map(convert_genbank, genbank_path_list))
    
    #print len(result)
    #print result.keys()
    
    handle = open(outfile, 'w')
    json.dump(result, handle)
    handle.close()
    #gc.collect()

def main():

    start = time.time()

    '''parser = argparse.ArgumentParser(description='Download precomputed HGT files or parse existing ones into a python dictionary that is stored in JSON.')

    parser.add_argument("-f", "--folder", dest="folder", metavar="FOLDER", default='./',
                help="Folder where results will be downloaded or read from, depending on the mode selected. Default is the folder where the script is run from.")
                
    parser.add_argument("-o", "--outfile", dest="outfile", default='hgt.json', metavar="FILE",
                help="Location where the program will store the json formatted result of the program. If omitted the program will use the folder plus 'HGT.json'.")


    parsed_args = parser.parse_args()
    folder = parsed_args.folder
    outfile = folder + parsed_args.outfile'''
    
    
    return_genbank_paths_of_interest('/home/dave/Desktop/all_genbank', './genbank_pathway_lists/', './genbank_pathway_lists/nc_filter_file.txt')
    
    #genbank_path_list = [i.strip() for i in open('./genbank_pathway_lists/filtered_list_genbank_paths.txt').readlines()]
    genbank_path_list  = [i.strip() for i in open('./genbank_pathway_lists/all_genbank_paths.txt').readlines()]
    #for fname in genbank_path_list:
    #    convert_genbank(fname)
    parallel_convert_genbank(genbank_path_list)

    
    #handle = open(outfile, 'w')
    #json.dump(hgt_dict, handle)
    #handle.close()


    print time.time() - start


if __name__ == '__main__':
    main()
