#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

# Copyright(C) 2013 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment

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
    
    # if the directory that the user specifies does not exist, then the program makes it for them.
    folder = parsed_args.folder  
    if not os.path.isdir(parsed_args.folder):
        os.makedirs(parsed_args.folder)
        
    outfile = parsed_args.outfile
    do_protein = parsed_args.protein
        
    return infile, outfile, folder, num_proc, do_protein
        

#################################################################
#### so yeah.... not gonna do this right now...             #####
#### I feel there is no reason to, but keeping              #####
#### the code in here just in case i am wrong about that    #####
#################################################################

# I need to do an analysis on gc content, and gc skew.
# currently this will return the  organism gc, (mean, SD, varience) of the GC in coding regions
# the results will we retained in a file the calling program opens tab delimited

def GCAnalysis(NC, organism, gc_list, seq, outfile):
    handle = open(outfile, 'a')
    organism_gc = "%3.2f" % GC(seq)
    mean = "%3.2f" % numpy.mean(gc_list)
    #mode = "%5.2f" % numpy.mode(gc_list)
    var = "%5.2f" % numpy.var(gc_list)
    std = "%5.2f" % numpy.std(gc_list)
    #print organism_gc
    handle.write('\t'.join([NC, organism, organism_gc, mean, var, std]) + NEW_LINE)
    

# take the genbank file specified by genbank path, and save the customized result file in the db_directory folder
#def convert_genbank(genbank_path, db_directory, error_fname): #, gc_outfile = 'gc_analysis.txt'):
def convert_genbank(genbank_tuple):
    genbank_path, db_directory, error_fname, do_protein = genbank_tuple
    record_list = []
    seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
    accession = seq_record.annotations['accessions'][0]
    organism = seq_record.annotations['organism'].replace(' ', '_')
    print 'organism ', organism
    err_log = []
    gc_list = [] # no need for this right now, but leaving in
    # loop over the genbank file
    for fnum, feature in enumerate(seq_record.features):
        err_flag = False
        if feature.type == 'CDS':
	    start = feature.location._start.position
            stop = feature.location._end.position
            strand = feature.strand
            seq = seq_record.seq[start:stop]
            gc = GC(seq)
            gc_list.append(gc)
            gc = "%3.2f" % gc
            
            if do_protein:
                seq = seq.translate()
                
            try:
                locus = feature.qualifiers['locus_tag'][0]
            except:
                try:
                    locus = feature.qualifiers['gene'][0]
                except:
                    locus = 'error'
                    print "Error in the organism %s with NC # %s" % (organism, accession)
                    err_flag = True
                    err_log.append([organism, accession])
            if len(seq) == 0:
                pass
            else:
                if 'gene' in feature.qualifiers:
                    gene = feature.qualifiers['gene'][0]
                    record_list.append(SeqRecord(seq, id = '|'.join([accession, organism, locus, gene, str(start), str(stop), str(strand), gc]).replace(' ', ''),
                       description = ''))
                else:
                    record_list.append( 
                      SeqRecord(seq, id = '|'.join([accession, organism, locus, 'unknown', str(start), str(stop), str(strand), gc]).replace(' ', ''),
                       description = ''))
    #if os.path.isfile(gc_outfile):
    #    os.remove(gc_outfile)
    #GCAnalysis(accession, organism, gc_list, seq_record.seq, gc_outfile)    
    handle = open(error_fname, 'a')
    for i in err_log:
        handle.write('\t'.join(i) + '\n')
        handle.close()
    if not err_flag:
        outpath = db_directory + os.path.splitext(os.path.basename(genbank_path))[0] + '.ffc'
        print outpath
        out_handle = open(outpath,"w")
        SeqIO.write(record_list, out_handle, "fasta")
        out_handle.close()
        
    if do_protein:
        cmd = "formatdb -i %s -p T -o F" % (outpath)
    else:    
       cmd = "formatdb -i %s -p F -o F" % (outpath)
    os.system(cmd)
        
        
    return outpath, err_flag


def parallel_convert_genbank(flist, folder, outfile, num_proc, error_fname, do_protein):
    pool = Pool(processes = num_proc)
    #result = dict(pool.map(convert_genbank, genbank_path_list))
    
    #convert_genbank(flist[0], folder, error_fname)
    
    tuple_list = [(i, folder, error_fname, do_protein) for i in flist]
    
    result = dict(pool.map(convert_genbank, tuple_list))
    
    handle = open(outfile, 'w')
    tmp = [i for i in result.keys() if not result[i]]
    handle.write('\n'.join(tmp))
    handle.close()


def main():
    
    start = time.time()
    
    parser = argparse.ArgumentParser(description='Convert a list of genbank files (full pathway) into a BLAST searchable database. Will save a list of created databases as a list of their pathways.')
                
    parser.add_argument("-i", "--infile", dest="infile", default='', metavar="FILE", required=True,
                help="A file that contains the full path to the every genbank file to be evaluated.")
                
    parser.add_argument("-o", "--outfile", dest="outfile", default='./blast_database_list.txt', metavar="FILE",
                help="Output file that contains the pathway to each of the newly created databases.")
    
    parser.add_argument("-f", "--folder", dest="folder", metavar="FOLDER", default='./db/',
                help="Folder where the databases that are created using formatdb will be stored. Default is the folder where the script is run from, but I do not reccomend using this option.")
    
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
    
    parser.add_argument("-p", "--protein", dest="protein", default=False, action='store_true',
                help="Folder where the databases that are created using formatdb will be stored. Default is the folder where the script is run from, but I do not reccomend using this option.")
    
                
    #parser.add_argument("-e", "--err_file", dest="err_file", default='./formatdb_err.log', metavar="FILE",
    #            help="Error file that contains the pathway to each of the newly created databases.")

    parsed_args = parser.parse_args()
    
    infile, outfile, folder, num_proc, do_protein = check_options(parsed_args)
    
    print infile, outfile, folder, num_proc, do_protein
        
    flist = [i.strip() for i in open(parsed_args.infile).readlines()] 
    
    # I might make this user selectable, but screw it for now.
    error_fname = "./error_log.txt"
    if os.path.isfile(error_fname):
        os.remove(error_fname)
    
    parallel_convert_genbank(flist, folder, outfile, num_proc, error_fname, do_protein)
    
    #convert_genbank(flist[0], folder, error_fname)
    
    print time.time() - start

    # A successful command could look like this:
    # ./format_db.py -i genbank_pathway_lists/filtered_list_genbank_paths.txt 
    
if __name__ == '__main__':
    main()
