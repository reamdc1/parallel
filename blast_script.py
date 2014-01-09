#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import argparse

# Copyright(C) 2013 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment

# This code right now only deals with protein, but I will add functionality later for nucleotides. 
# Just moving the project along here.
def do_parallel_blast(arg_tuple):
    db, query_file, blast_result_folder, num_processors, eval_threshold = arg_tuple
    out_file = "%s%s_prot.txt" % (blast_result_folder, db.split('/')[-1].split('.')[0])
    cmd = "blastall -p tblastn -a %i -i %s -d %s -e %s -o %s -m 9" % (num_processors, query_file, db, eval_threshold, out_file)
    os.system( cmd )

def check_options(parsed_args):
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = parsed_args.num_proc
        
    if os.path.exists(parsed_args.infile):
        infile = parsed_args.infile
    else:
        print "The file %s does not exist." % parsed_args.infile
        sys.exit()
        
    if os.path.exists(parsed_args.query):
        query = parsed_args.query
    else:
        print "The file %s does not exist." % parsed_args.query
        sys.exit()
    
    # if the directory that the user specifies does not exist, then the program makes it for them.
    folder = parsed_args.folder  
    if not os.path.isdir(parsed_args.folder):
        os.makedirs(parsed_args.folder)
    # should offer some sanity here for the input, but it is missing right now. call it a to do item    
    e_val = parsed_args.eval
    return infile, query, folder, num_proc, e_val


def parallel_blast(infile, query, folder, num_proc, e_val = '1e-10'):
    # you kinda have to trust me here, but having blast run on as many threads per CPU as you have total processors is fastest
    # I have no idea why this is... ugh.
    blast_arg_list = [(i.strip(), query, folder, num_proc, e_val) for i in open(infile).readlines()]
    pool = Pool(processes = num_proc)
    pool.map(do_parallel_blast, blast_arg_list)

def main():
    
    start = time.time()
    
    parser = argparse.ArgumentParser(description="Conduct a BLAST search over a list of  BLAST searchable databases using a common query file. The program will save the results in a folder designated by the user or the default './blast_result/'.")
                
    parser.add_argument("-i", "--infile", dest="infile", default='', metavar="FILE", required=True,
                help="A file that contains the path to every organism database that you are interested in.")
    
    parser.add_argument("-f", "--folder", dest="folder", metavar="FOLDER", default='./blast_result/',
                help="Folder where the BLAST results will be stored. Default is the folder './blast_result/'.")
    
    parser.add_argument("-q", "--query", dest="query", default='./blast_database_list.txt', metavar="FILE",
                help="A file that contains the BLAST query for every gene of interest in the dataset.")
    
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")

    parser.add_argument("-e", "--eval", dest="eval", default='1e-10', metavar="FLOAT",
                help="eval for the BLAST search.")
    
    parsed_args = parser.parse_args()
    
    infile, query, folder, num_proc, e_val = check_options(parsed_args)
    
    parallel_blast(infile, query, folder, num_proc, e_val)

    print time.time() - start

    # ./blast_script.py -i blast_database_list.txt -f ./blast_result/ -q ./fasta_generation/protein_matches.fa 
    
if __name__ == '__main__':
    main()

