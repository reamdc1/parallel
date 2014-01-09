from multiprocessing import Pool
import numpy
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
#from homolog import *
from homolog2 import *
from collections import defaultdict
import simplejson as json

# Global variables... blah blah, fast and easy rules the day here
NUM_PROCESSORS = os.sysconf("SC_NPROCESSORS_CONF")


# Ok so this bit of code is designed to find all files in genbank that are stored locally and report their entire
# path, saved in ./AllGenbankPath.txt. Next it reads in the NC_number if the genomes we are interested in. It matches 
# this to the file path, then extracts the home folder of said genome. The genome folder is then copied into a small
# directory. The copying of files should be removed if we are planning on working on the whole data set.
# This function additionally downloads IslandPath HGT information, parses the files into a dictionary, 
# then saves them to a JSON file for fast recovery and portability.
def startup_crap(copy_folders = True, download_hgt = True):
    res = return_recursive_dir_files('/home/dave/Desktop/all_genbank')
    #res = return_recursive_dir_files('./proteobacteria_genbank')
    #handle = open('./AllGenbankPath.txt', 'w')
    handle = open('./recover/AllGenbankPath2.txt', 'w')
    handle.write('\n'.join(res))
    handle.close()
    
    NC_list = [i.strip() for i in open('./phylo_order.txt').readlines()]
    
    result_list = [i for i in res if i.split('/')[len(i.split('/')) -1].split('.')[0] in NC_list]
    
    for item in result_list:
        folder = os.path.dirname(item)
        #print folder
        if copy_folders:
            os.system( "cp -R %s/ ./proteobacteria_genbank/" % folder)
            
    if download_hgt:
        #download_fiona_hgt('./hgt/')        
        os.system("./hgt.py -d -f ./hgt/ -o hgt.json")
    




def download_remote_file(url, outfile, block = 4096): # low memory implementation
    """Overview: This function downloads the remote file specified in url and saves its as outfile
       while limiting the memory footprint for large files.
       
       Return: a flie downloaded from the url location and saved as outfile
       
       Deafult: block size of 4096, which is convenient to keep memory usage low 
    """
    import urllib2
    dl_file = urllib2.urlopen(url)
    handle = open(outfile, 'w')
    while 1:
        data = dl_file.read(block)
        if data:
            handle.write(data)
        else:
            break
    handle.close()

# This function will download the HGT files from Fiona Brinkman. This is just a quick thing to update the 2files if I need this
# to run later.

def download_fiona_hgt(result_folder): # i could make this parallel and i might.
    
    # This file specifically does not tell me what method was used to determine the data, so it is not very informative for the purpose we want.
    #download_remote_file('http://www.pathogenomics.sfu.ca/islandviewer/download_all.php?tool=islandviewer&type=csv', './hgt/all_gis_islandviewer.csv')
    
    download_remote_file('http://www.pathogenomics.sfu.ca/islandviewer/download_all.php?tool=islandpick&type=csv', result_folder + 'all_gis_islandpick.csv')
    
    download_remote_file('http://www.pathogenomics.sfu.ca/islandviewer/download_all.php?tool=islandpath_dimob&type=csv', result_folder + 'all_gis_islandpath_dimob.csv')
    
    download_remote_file('http://www.pathogenomics.sfu.ca/islandviewer/download_all.php?tool=sigi_hmm&type=csv', result_folder + 'all_gis_sigi_hmm.csv')
    

def parse_hgt_files(flist):
    result = {}
    for fname in flist:
        for element in [i.strip() for i in open(fname).readlines()[1:]]:
            accession, start, stop, size, program = [i.split('.')[0] for i in element.split(',')]
            if accession in result.keys():
                result[accession].append((int(start), int(stop), int(size), program))
            else:
                result.update({accession:[(int(start), int(stop), int(size), program)]})

    return result

#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def return_recursive_dir_files(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


def create_distance_dict(dist_file):
    result = {}
    
    for line in [i.strip() for i in open(dist_file).readlines()]:
        NC, dist = line.split('\t')
        result.update({NC: float(dist)})
    return result


# this is not a piece of code that should remain. i have to include it for testing since i am developing this module seperatly
# from any main fork.  


def convert_genbank(genbank_tuple):
    genbank_path, hgt_dict, db_directory = genbank_tuple
    #print 'hooty hoo', genbank_path#, hgt_dict, db_directory
    record_list = []
    seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
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
    if re.search("hromosome", description) is not None: # we have a segmented chromosome (more than one chromosome)
        seq_type = 'chromosome'
        seq_name = description.split('hromosome')[1].split(',')[0].strip()
    elif re.search("lasmid", description) is not None: # we have a plasmid
        seq_type = 'plasmid'
        seq_name = description.split('lasmid')[1].split(',')[0].strip()
    else:
        seq_type = 'chromosome'
        seq_name = 'complete_genome'

    #print 'path', genbank_path, 'source folder', source_folder, 'organism ', organism
       
    # loop over the genbank file
    for fnum, feature in enumerate(seq_record.features):
        if feature.type == 'CDS':
            start = feature.location._start.position
            stop = feature.location._end.position
            strand = feature.strand
            #print 'keys = ', feature.qualifiers.keys()
            
            # This next block of code should determine HGT, and will only provide
            hgt = []

            if accession in hgt_dict.keys():
                hgt = ['NONE']
                hgt_list = hgt_dict[accession]
                #print len(hgt_list)
                for item in hgt_list:
                    hgt_start, hgt_stop, hgt_size, hgt_program = item
                    if start >= hgt_start and start <= hgt_stop:
                        #print NC, org, predicted, hgt_program, 'hgt suspected!'
                        if hgt[0] == 'NONE':
                            hgt = [hgt_program]
                        else:
                            hgt.append(hgt_program)
            else:
                hgt = ['NotEvaluated']
                print accession, 
            if hgt[0] != 'NONE' and hgt[0] != 'NotEvaluated':
                #print 'hgt ', ':'.join(sorted(list(set(hgt))))
                hgt = ':'.join(sorted(list(set(hgt))))
                print accession, hgt
            else:
                hgt = hgt[0]
            try:
                locus = feature.qualifiers['locus_tag'][0]
            except:
                locus = feature.qualifiers['gene'][0]
            seq = seq_record.seq[start:stop]
            gc = "%2.1f" % GC(seq)

            if len(seq) == 0:
                pass
            else:
                if 'gene' in feature.qualifiers:
                    gene = feature.qualifiers['gene'][0]
                    #record_list.append(SeqRecord(seq, id='ref|' + '|'.join([accession, locus, gene, str(start), str(stop), str(strand), organism, gc]),
                    #   description = ''))  
                    record_list.append(SeqRecord(seq, id='ref|' + '|'.join([accession, locus, gene, str(start), str(stop), str(strand), organism, gc, hgt]),
                       description = '')) 
                else:
                    #record_list.append(SeqRecord(seq, id='ref|' + '|'.join([accession, locus, 'unknown', str(start), str(stop), str(strand), organism, gc]),
                    #   description = ''))
                    record_list.append(SeqRecord(seq, id='ref|' + '|'.join([accession, locus, 'unknown', str(start), str(stop), str(strand), organism, gc, hgt]),
                       description = ''))

    outpath = os.path.splitext(os.path.basename(genbank_path))[0] + ".ffc"
    SeqIO.write(record_list, open(db_directory + outpath,"w"), "fasta")
    
    cmd = "formatdb -i %s -p F -o F" % (db_directory + outpath)
    os.system(cmd)

DB_FOLDER = './recover/db_folder/'
def make_bast_db(genbank_path_list_file, num_proc):

    hgt_dict = parse_hgt_files(return_recursive_dir_files('./hgt'))
    #print len(hgt_dict.keys())
    
    flist = [(i.strip(), hgt_dict, DB_FOLDER) for i in open(genbank_path_list_file).readlines()]
    pool = Pool(processes = num_proc)
    result = pool.map(convert_genbank, flist)

# The name sucks, forgive me. This function has the purpose of taking each organism and its phylogenetic distance (from E. coli)
# and returning a dictionary keyed to that organism. the dictionary elements will be a list that contains the in-order closest
# relatives.  
def return_closest_relatives(distance_file):
    dist_dict = create_distance_dict(distance_file)
    nc_list = [i.split('\t')[0] for i in open(distance_file).readlines()]
    #print len(sorted(itertools.combinations(nc_list, 2)))
    sweet_list = sorted(itertools.permutations(nc_list, 2)) # stupid name, sorry
    
    result = {}
    
    for item in sweet_list:
        primary_org, compare_org = item
        p_dist = dist_dict[primary_org]
        c_dist = dist_dict[compare_org]
        dist = math.fabs(p_dist-c_dist)
        if primary_org not in result.keys():
            result.update({primary_org: [(compare_org, dist)]})
        else:
            result[primary_org].append((compare_org, dist))
            
            
    for item in result.keys():
        sort_tmp = sorted(result[item], key=lambda t: t[1])
        #print sort_tmp
        result.update({item: sort_tmp})

    return result
        
    
# This function creates a dictionary indexed by locus from the input genbank file
# and for my purposes now, it will index genes based on their annotation in genbank
# seq_type will allow us to determine the type of sequence returned in for the dict. default
# will be amino acid because this is a lot less noisy.
def return_genbank_dict(gb_file, key = 'annotation', seq_type = 'amino_acid'):
    """Overview: This function will return a dictionary generated from a genbank file with key value supplied by caller.
       Returns: A dictionary created by the supplied genbank file (gb_file) indexed off the key value supplied.
       Default: The deafult key is locus, and this is generally the most useful key type since it is garanteed to be 
       unique within the genbank file. This condition is not necessarily true for any other attribute.
   """
    result = {}
    seq_record = SeqIO.parse(open(gb_file), "genbank").next()
    accession = seq_record.annotations['accessions'][0].split('.')[0]
    common_name = seq_record.annotations['organism'].replace(' ', '_')
    result.update({'accession': accession})
    result.update({'common_name': common_name})
    cnt = 0
    # loop over the genbank file
    unk_cnt = 1
    for fnum, feature in enumerate(seq_record.features):
        # here i simply check the gene coding type, and identify them in a way that can be used later.
        if feature.type == 'CDS' or feature.type == 'ncRNA' or feature.type == 'tRNA' or feature.type == 'mRNA' or feature.type == 'rRNA':
	    start = feature.location._start.position
            stop = feature.location._end.position
            strand = feature.strand
            synonyms = 'NONE'
            if 'gene_synonym' in feature.qualifiers:
                synonyms = ':'.join(feature.qualifiers['gene_synonym'][0].replace(' ', '').split(';'))
            try:
                locus = feature.qualifiers['locus_tag'][0]
            except:
                try:
                    locus = feature.qualifiers['gene'][0]
                except:
                    locus = ''
                    print 'No locus associated. This should never be invoked meaning you are proper fracked. (The gbk file has an error).'
            try: 
                gene = feature.qualifiers['gene'][0]
            except:
                gene = 'unknown'
            try:
                seq = feature.qualifiers['translation']
                seq_type = 'Protein'
            except:
                cnt = cnt + 1
                seq = seq_record.seq[start:stop]
                seq_type = feature.type
                if feature.type == 'CDS':
                    seq_type = 'Pseudo_Gene'
            gc = "%2.1f" % GC(seq_record.seq[start:stop])
            #method = "exact"
            if key == 'locus':
                result.update({locus: (locus, gene, seq, seq_type, synonyms)})
            elif key == 'annotation':
                if gene == 'unknown':
                    new_gene = 'unknown_' + str(unk_cnt)

                    header = ">" + '|'.join([accession, common_name, locus, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
                    result.update({new_gene: [header, ''.join(seq)]})
                    unk_cnt +=1
                else:
                    header = ">" + '|'.join([accession, common_name, locus, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
                    result.update({gene: [header, ''.join(seq)]})

    #print 'The number of non-protein regions in %s is: %i.' % (common_name, cnt)
    return result


# This function will allow me to do the main work of make_operon_fasta, but allow parallel
# processing. I will have to make a parallel array, which will take some time to learn.  i 
# will keep this stub for later to implement. 
def parallel_operon_fasta(genome):
    organism = genome.split('/')[-1].split('.')[0]
    organism_dict_for_recovery = {}
    org_dict = return_genbank_dict(genome)
    organism_dict_for_recovery.update({organism: org_dict})
    return (organism, org_dict)
        

# This function will make a BLAST query from the parsed operon file.  The default behavior of the function is 
# to make a query file from all genes in all organisms that are annotated. Later the results will be sorted based 
# the needs of the programmer. The defaulted variables allow a single operon to be chosen individually.  The program
# will also store teh results of this function in a folder titled blast_query_files, in the recovery folder.
def make_operon_fasta(operon_file, result_folder = './recover/blast_query_files/', operon_considered = 'ALL', filter_file = 'NONE'):

    #print "operon_file", operon_file
    distance_file = './recover/Dist_test.txt'
    dist_org_list = [i.split('\t')[0] for i in open(distance_file).readlines()]
    print dist_org_list
    
    res = return_recursive_dir_files('./proteobacteria_genbank')
    genome_of_interest_list = [i for i in res if i.split('/')[len(i.split('/')) -1].split('.')[0] in dist_org_list]
    #print len(genome_of_interest_list), genome_of_interest_list
    
    organism_dict_for_recovery = {}
    
    if filter_file == 'NONE': # look for everything that we have a distance for. This is goign to be the best option
        pool = Pool(processes = NUM_PROCESSORS)
        organism_dict_for_recovery = dict(pool.map(parallel_operon_fasta, genome_of_interest_list))

    else: # TO DO: put  code in here to parse the filter file and then enforce that an NC number be in this file before evaluation
        pass
        
    exact_match = []
    rna_exact_match = []
    missing_list = []

    for operon_line in [i.strip() for i in open(operon_file).readlines()]:
        operon = operon_line.split('\t')[0]
        if operon_considered == 'ALL':
            gene_list = operon_line.split('\t')[1:]
            for org in organism_dict_for_recovery.keys():
                #print "got here"
                for gene in gene_list:
                    if gene in organism_dict_for_recovery[org].keys():
                        #print organism_dict_for_recovery[org][gene]
                        if organism_dict_for_recovery[org][gene][0].split('|')[7] == 'Protein':
                            
                            exact_match.append(organism_dict_for_recovery[org][gene][0])
                            exact_match.append(organism_dict_for_recovery[org][gene][1])
                            #print organism_dict_for_recovery[org][gene][0]
                        else:
                            print organism_dict_for_recovery[org][gene][0]
                            #exact_match.append(organism_dict_for_recovery[org][gene][1])
                        #print organism_dict_for_recovery[org][gene][1]
                    else: # The gene is missing, we will look for it at a later stage in the program
                        item = '\t'.join([org, gene])
                        missing_list.append(item)
            #print "all operons will be used (but this is not working yet)."
        elif operon_considered == operon:
            gene_list = operon_line.split('\t')[1:]
            for org in organism_dict_for_recovery.keys():
                print org
                for gene in gene_list:
                    if gene in organism_dict_for_recovery[org].keys():
                        #print organism_dict_for_recovery[org][gene]
                        exact_match.append(organism_dict_for_recovery[org][gene][0])
                        exact_match.append(organism_dict_for_recovery[org][gene][1])
                    else: # The gene is missing, we will look for it at a later stage in the program
                        item = '\t'.join([org, gene])
                        missing_list.append(item)
        handle = open(result_folder + 'exact_matches.fa', 'w')
        handle.write('\n'.join(exact_match))
        handle.close()
        handle = open(result_folder + 'missing_operon_genes.txt', 'w')
        handle.write('\n'.join(missing_list))
        handle.close()

        #else:
        #    pass

def parallel_blast(db_folder, query_file, blast_result_folder, num_processors, eval_threshold):
    tmp = return_recursive_dir_files(db_folder)
    db_list = []
    for item in tmp:
        fname, file_ext = os.path.splitext(item)
        if file_ext == '.ffc':
            db_list.append((item, query_file, blast_result_folder, num_processors, eval_threshold))
        else:
            pass
    pool = Pool(processes = NUM_PROCESSORS)
    pool.map(do_parallel_blast, db_list)

def do_parallel_blast(arg_tuple):
    db, query_file, blast_result_folder, num_processors, eval_threshold = arg_tuple
    out_file = "%s%s_prot.txt" % (blast_result_folder, db.split('/')[len(db.split('/')) - 1].split('.')[0])
    cmd = "blastall -p tblastn -a %i -i %s -d %s -e %s -o %s -m 9" % (num_processors, query_file, db, eval_threshold, out_file)
    os.system( cmd )

# If we do not wish to filter any of the results, nc_filter_file should be set to the value 'NONE'
# I am only workign with the _prot right now, but we will have to deal with the rna thing too, later
def parallel_blast_parse(nc_filter_file, blast_folder, result_folder, operon_file, num_procs = NUM_PROCESSORS):
    # this block from here to the first pool.map line is unnecessary. I am using this as a function tracing, but not for analysis. this code block should be removed
    # since it serves no real purpose in the long term.
    if nc_filter_file == 'NONE':
        filtered_files = [(i, './recover/blast_parse/raw_info/')  for i in return_recursive_dir_files(blast_folder)]
    else:
        nc_filter = [i.strip() for i in open(nc_filter_file).readlines()]
        filtered_files =  [(i, './recover/blast_parse/raw_info/')  for i in return_recursive_dir_files(blast_folder) if i.split('/')[len(i.split('/'))-1].split('_p')[0] in nc_filter]

    pool = Pool(processes = NUM_PROCESSORS)
    pool.map(initial_parse_blast, filtered_files)
    
    if nc_filter_file == 'NONE':
        filtered_files2 = [(i, result_folder) for i in return_recursive_dir_files(blast_folder)]
    else:
        nc_filter = [i.strip() for i in open(nc_filter_file).readlines()]
        filtered_files2 =  [(i, result_folder)  for i in return_recursive_dir_files(blast_folder) if i.split('/')[len(i.split('/'))-1].split('_p')[0] in nc_filter]

    pool = Pool(processes = num_procs)
    pool.map(parse_blast, filtered_files2)
    
    
# this is to be used with parallel_blast_parse ().  
# this is to make the first set of debugging files, sorry about it, but just have to deal with this while i
# am debugging things.
def initial_parse_blast(parse_tuple):
    fname, result_folder = parse_tuple
    hit_list = [i.strip() for i in open(fname).readlines() if i[0] != '#']

    #outfile = './recover/blast_parse/raw_info/' + fname.split('/')[len(fname.split('/'))-1].split('_p')[0] + '.txt'
    outfile = './recover/blast_parse/raw_info/' + fname.split('/')[-1].split('_p')[0] + '.txt'


    handle = open(outfile, 'w')
    handle.write('\n'.join(hit_list))
    handle.close()


# this function does the heavy lifting for BLAST parsing.     
def parse_blast(parse_tuple):
    fname, result_folder = parse_tuple
    
    # this will be locus keyed raw data.... at least for now.....
    locus_hit_dict = {}
    homolog_locus_hit_dict = {}
    
    hit_list = [i.strip() for i in open(fname).readlines() if i[0] != '#']
    
    for item in hit_list:
        locus = item.split('\t')[1].split('|')[2]
        if locus in locus_hit_dict.keys():
            locus_hit_dict[locus].append(item)
        else:
            locus_hit_dict.update({locus: [item]})    
        
        if locus not in homolog_locus_hit_dict.keys():
            homolog_locus_hit_dict.update({locus: []})
        #homolog_locus_hit_dict[locus].append(res)
        homolog_locus_hit_dict[locus].append(blast_hit_to_homolog(item))
    #print "homolog_locus_hit_dict", len(homolog_locus_hit_dict)

    
    result = []
    for locus in locus_hit_dict.keys():
        result.append(locus)
        for line in locus_hit_dict[locus]:
            result.append(line)
            
 
    #outfile = result_folder + fname.split('/')[len(fname.split('/'))-1].split('_p')[0] + '.txt'
    #outfile1 = result_folder + fname.split('/')[len(fname.split('/'))-1].split('_p')[0] + '_homolog1.txt'
    outfile1 = result_folder + fname.split('/')[len(fname.split('/'))-1].split('_p')[0] + '.txt'
    
    #handle = open(outfile, 'w')
    #handle.write('\n'.join(result))
    #handle.close()

    homolog_result = []
    
    #for locus in homolog_locus_hit_dict.keys():
    #    homolog_result.append(locus)
    #    for hlog in homolog_locus_hit_dict[locus]:
    #        homolog_result.append(hlog.ret_str())
    
    tmp = filter_homolog_locus_dict(homolog_locus_hit_dict)
    
    self_homolog_dict = json.load(open('./recover/self_homologs_json.txt'))
    tmp2 = better_homolog_locus_filter(homolog_locus_hit_dict, self_homolog_dict)
    print tmp2
    for locus in sorted(tmp.keys()):
        #homolog_result.append(locus)
        homolog_result.append(tmp[locus].ret_str())

    handle1 = open(outfile1, 'w')
    handle1.write('\n'.join(homolog_result))
    handle1.close()


# this function will take the list of homologs obtained per locus and return a dict that has been filtered per locus
# in an intelligent manner. the goal is to have the best represenative BLAST hit for a given loci.


# this is something that i am goign to have to re-write.  i have to allow multiple hits per locus, despite my intenet not.
def filter_homolog_locus_dict(hlog_dict):
    result = {}

    for locus in sorted(hlog_dict.keys()):
        h_list = hlog_dict[locus]
        best_homolog = h_list[0]
        for h in h_list[1:]: # the first in the list is already the default best, yahoo.
            if h.method() == 'exact':
                best_homolog = h
                break
            elif h.percent_ident > best_homolog.percent_ident:
                best_homolog = h

        result.update({locus: best_homolog})
        
    return result
            

def better_homolog_locus_filter(hlog_dict, self_hlog_dict):

    for locus in sorted(hlog_dict.keys()):
        h_list = hlog_dict[locus]
        blast_prediction_count = defaultdict(int)
        #print locus, " Length : " , len(h_list)
        for item in h_list:
            blast_prediction_count[ item.predicted_gene()] += 1 #, item.predicted_gene()
            #item.Print()
        #for item in h_list:
        #    item.Print()
        print blast_prediction_count.viewitems()

    return 1


# this function will take the list of homologs obtained per locus and return a dict that has been filtered per locus
# in an intelligent manner. the goal is to have the best represenative BLAST hit for a given loci.
def filter_homolog_locus_dict_old(hlog_dict):
    result = {}
    for locus in sorted(hlog_dict.keys()):
        h_list = hlog_dict[locus]
        best_homolog = h_list[0]
        for h in h_list[1:]: # the first in the list is already the default best, yahoo.
            if h.method() == 'exact':
                best_homolog = h
                break
            elif h.percent_ident > best_homolog.percent_ident:
                best_homolog = h

        result.update({locus: best_homolog})
        
    return result



# This function is a parser for converting a BLAST hit into a Homolog class object. 



# To clarify what is going on here: BLAST hits run from the -m 9 tabular with comment lines.  I filter out the comment lines, but makes everything 
# more readable in my opinion.  

# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score


def blast_hit_to_homolog(hit):
    field = hit.split('\t')
    query = field[0].split('|')
    source_accession = query[0]
    source_common = query[1]
    source_locus = query[2]
    predicted_gene = query[3]
    source_start = query[4]
    source_start = query[5]
    source_strand = query[6]
    product_type = query[7]
    synonyms = query[8]
    synonym_list = synonyms.split(':')
    source_gc = query[9]
    accession, locus, gene, start, stop, strand, organism, gc, hgt = field[1].split('|')[1:]
    percent_ident = field[2]
    
    alignment_length = field[3]
    mismatches = field[4]
    gap_openings = field[5]
    query_start = field[6]
    query_end = field[7]
    subject_start = field[8]
    subject_end = field[9]
    e_val = field[10]
    bit_score = field[11]
    alignment_length = field[3]
    #operon_name = gene_to_operon_dict[predicted_gene]
    method = 'BLAST' # might have to change this to the version of blast that is used... not sure there.
    if source_accession == accession and source_locus == locus: # eliminates problem of homologous genes within a single operon
        method = "exact"
    
        
    return Homolog(accession, organism, locus, gene, predicted_gene, synonyms, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type, alignment_length, method, source_accession, source_common, source_locus, source_start, hgt)


def make_gene_to_operon_dict(operon_file, parsed = True):
    result = {}
    if parsed:
        lst = [i.strip() for i in open(operon_file).readlines()]
    else: #This is where i need to call the function that parses the regulonDB file/downloads it.... 
        pass

    for item in lst:
        tmp = item.split('\t')
        operon = tmp[0]
        for gene in tmp[1:]:
            result.update({gene:operon})
    return result
    

def operon_sort(operon_file, blast_parse_folder, nc_filter_file = 'NONE'):
    result = {}
    gene_to_operon_dict = make_gene_to_operon_dict(operon_file)
    for operon in [i.split('\t')[0] for i in open(operon_file).readlines()]:
        result.update({operon:[]})
    if nc_filter_file == 'NONE':
        filtered_files = [i for i in return_recursive_dir_files(blast_parse_folder)]
    else:
        nc_filter = [i.strip() for i in open(nc_filter_file).readlines()]
        # might have to change the '_p' bit below... 
        filtered_files =  [i for i in return_recursive_dir_files(blast_folder) if i.split('/')[len(i.split('/'))-1].split('_p')[0] in nc_filter]
        
    for fname in filtered_files:
        for line in open(fname).readlines():
            res = Homolog.from_file(line)
            operon = gene_to_operon_dict[res.predicted_gene()]
            result[operon].append(res)
    return result



# this function will return a dictionary of operon keyed off the operon name with data values in the form
# of a list of homologs which are homologous. ex. [abcA, abcB]

# i have conspiciously removed RNA stuff here.  that is on purpose, it was causing errors, and i don't have time to debug it now.
def return_self_homolog_dict(operon_list = './recover/operon_name_and_genes.txt', prot_file = './recover/operon_protein_query.fa', rna_file = './recover/operon_rna_query.fa', outfile = './recover/self_homologs_json.txt'):
    # makes a dictionary keyed by operon name and a list of the gene contained by the operon
    operon_dict = {}
    for line in [i.strip() for i in open(operon_list).readlines()]:
        tmp = line.split('\t')
        operon_dict.update({tmp[0]:tmp[1:]})
        
    # set up databases for the different types of genes
    # for proteins -p must be set to true    
    cmd = "formatdb -i %s -p T -o F" % (prot_file)
    os.system(cmd)
    # for RNA genes -p must be set to false 
    #cmd = "formatdb -i %s -p F -o F" % (rna_file)
    #os.system(cmd)

    
    # blast each set of genes against itself
    cmd = "blastall -p blastp -a %i -i %s -d %s -e %s -o %s -m 9" % (os.sysconf("SC_NPROCESSORS_ONLN"), prot_file, prot_file, '1e-10', './recover/self_prot.txt')
    os.system( cmd )
    #cmd = "blastall -p blastn -a %i -i %s -d %s -e %s -o %s -m 9" % (os.sysconf("SC_NPROCESSORS_ONLN"), rna_file, rna_file, '1e-10', './recover/self_rna.txt')
    #os.system( cmd )

    
    # in this next section i will read in the resulting blast results, and construct a dictionary which will be keyed off gene name and provide a list
    # of homologs from the operon set. This list will allow the program to filter out spurious results. We will miss fusions of homologous genes, but
    # hopefully this will be a rare event in our dataset, untill this can be revised
    
    lst = [i.strip() for i in open('./recover/self_prot.txt').readlines() if i[0] != '#']
    for line in [i.strip() for i in open('./recover/self_rna.txt').readlines() if i[0] != '#']:
        lst.append(line)
        
    result = {} 
    
    for line in lst:
        source, hit = line.split('\t')[0:2]
        source_annotation = source.split('|')[2]
        hit_annotation = hit.split('|')[2]
        # we have two genes in the test set that are homologous
        if source_annotation != hit_annotation:
            if source_annotation not in result.keys():
                result.update({source_annotation: [hit_annotation]})
            else:
                result[source_annotation].append(hit_annotation)
    print result  
    json_handle = open(outfile, 'w')
    json.dump(result, json_handle)
    json_handle.close()
    return result



def main():

    # This file contains the phylogenetic distances from E. coli, as determined from rpoD. (unless it is rpoB). 
    # Either way, it contains two values. NC number and distance from E. coli, tab delineated.
    distance_file = './recover/Dist_test.txt'

    start = time.time()
    
    # this function should be called once before you start running the rest. not every time. ugh.
    #startup_crap(True)
    startup_crap(False, False)

    #hgt_dict = parse_hgt_files(return_recursive_dir_files('./hgt'))
    #print len(hgt_dict.keys())
    
    ##make_bast_db('./recover/AllGenbankPath.txt', NUM_PROCESSORS)
    
    ##closest_relatives_dict = return_closest_relatives(distance_file)

    #make_operon_fasta('./recover/operon_name_and_genes.txt', './recover/blast_query_files/', 'atpIBEFHAGDC')
    make_operon_fasta('./recover/operon_name_and_genes.txt', './recover/blast_query_files/')
     
    eval_threshold = .001
    
    #parallel_blast('./recover/db_folder/', './recover/blast_query_files/exact_matches.fa', './recover/blast_result/', NUM_PROCESSORS, eval_threshold)
    
    #parallel_blast_parse('phylo_order.txt', './recover/blast_result/', './recover/blast_parse/intermediate/', './recover/operon_name_and_genes.txt', NUM_PROCESSORS)
    
    self_homolog_dict = return_self_homolog_dict()
    
    parallel_blast_parse('phylo_order.txt', './recover/blast_result/', './recover/blast_parse/intermediate/', './recover/operon_name_and_genes.txt', 1)
    
    #parallel_blast_parse('NONE', './recover/blast_result/', './recover/blast_parse/intermediate/', './recover/operon_name_and_genes.txt', NUM_PROCESSORS)
    
    # In this part, we now take the each blast_parsed file and make a huge frickin' result based on the name of an operon.  yahoo for school!
    
    ##operon_sort('./recover/operon_name_and_genes.txt', './recover/blast_parse/intermediate/')
    

    end = time.time()
    print end - start




if __name__ == '__main__':
    main()

