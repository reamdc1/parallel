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
from homolog import *

# Global variables... blah blah, fast and easy rules the day here
NUM_PROCESSORS = os.sysconf("SC_NPROCESSORS_CONF")


# Ok so this bit of code is designed to find all files in genbank that are stored locally and report their entire
# path, saved in ./AllGenbankPath.txt. Next it reads in the NC_number if the genomes we are interested in. It matches 
# this to the file path, then extracts the home folder of said genome. The genome folder is then copied into a small
# directory. The copying of files should be removed if we are planning on working on the whole data set.
def startup_crap(copy_folders = True):
    #res = return_recursive_dir_files('/home/dave/Desktop/all_genbank')
    res = return_recursive_dir_files('./proteobacteria_genbank')
    #handle = open('./AllGenbankPath.txt', 'w')
    handle = open('./recover/AllGenbankPath.txt', 'w')
    handle.write('\n'.join(res))
    handle.close()
    
    NC_list = [i.strip() for i in open('./phylo_order.txt').readlines()]
    
    result_list = [i for i in res if i.split('/')[len(i.split('/')) -1].split('.')[0] in NC_list]
    #print len(result_list)
    
    for item in result_list:
        folder = os.path.dirname(item)
        #print folder
        if copy_folders:
            os.system( "cp -R %s/ ./proteobacteria_genbank/" % folder)

    

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
DB_FOLDER = './recover/db_folder/'

def convert_genbank(genbank_path, db_directory  = DB_FOLDER):
   
    record_list = []
    seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
    nc_number = seq_record.name
    
    # genbank information about the sequence found in the file. 
    accession = seq_record.id.split('.')[0]
    print "got here", accession
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

    print 'path', genbank_path, 'source folder', source_folder, 'organism ', organism
    
    
    
    # loop over the genbank file
    for fnum, feature in enumerate(seq_record.features):
        if feature.type == 'CDS':
            start = feature.location._start.position
            stop = feature.location._end.position
            strand = feature.strand
            #print 'keys = ', feature.qualifiers.keys()
            try:
                locus = feature.qualifiers['locus_tag'][0]
            except:
                locus = feature.qualifiers['gene'][0]
            seq = seq_record.seq[start:stop]
            if len(seq) == 0:
                pass
            else:
                if 'gene' in feature.qualifiers:
                    gene = feature.qualifiers['gene'][0]
                    record_list.append(SeqRecord(seq, id="ref|%s|%s|%s|%s|%s|%s|%s" % (accession, locus, gene, str(start), str(stop), str(strand), organism),
                       description = ''))   
                else:
                    record_list.append( 
                      SeqRecord(seq, id="ref|%s|%s|%s|%s|%s|%s|%s" % (accession, locus, 'unknown', str(start), str(stop), str(strand), organism),
                       description = ''))

    outpath = os.path.splitext(os.path.basename(genbank_path))[0] + ".ffc"
    SeqIO.write(record_list, open(db_directory + outpath,"w"), "fasta")
    
    cmd = "formatdb -i %s -p F -o F" % (db_directory + outpath)
    os.system(cmd)


def make_bast_db(genbank_path_list_file, num_proc):
    flist = [i.strip() for i in open(genbank_path_list_file).readlines()]
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
    #print result['NC_000913']
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
                    
                    #result.update({new_gene: (accession, locus, gene, str(start), str(stop), str(strand), common_name, gc)})
                    header = ">" + '|'.join([accession, common_name, locus, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
                   
                    #result.update({new_gene: (accession, common_name, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc, method)})
                    result.update({new_gene: [header, ''.join(seq)]})
                    unk_cnt +=1
                else:
                    #result.update({gene: (accession, locus, gene, str(start), str(stop), str(strand), common_name, gc)})
                    #result.update({gene: (accession, common_name, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc, method)})
                    header = ">" + '|'.join([accession, common_name, locus, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
                    result.update({gene: [header, ''.join(seq)]})

    #print 'The number of non-protein regions in %s is: %i.' % (common_name, cnt)
    return result


# This function will allow me to do the main work of make_operon_fasta, but allow parallel
# processing. I will have to make a parallel array, which will take some time to learn.  i 
# will keep this stub for later to implement. 
def parallel_operon_fasta(genome):
    organism = genome.split('/')[len(genome.split('/')) -1].split('.')[0]
    print organism
    organism_dict_for_recovery = {}
    org_dict = return_genbank_dict(genome)
    organism_dict_for_recovery.update({organism: org_dict})
    #return organism_dict_for_recovery
    return (organism, org_dict)
        

# This function will make a BLAST query from the parsed operon file.  The default behavior of the function is 
# to make a query file from all genes in all organisms that are annotated. Later the results will be sorted based 
# the needs of the programmer. The defaulted variables allow a single operon to be chosen individually.  The program
# will also store teh results of this function in a folder titled blast_query_files, in the recovery folder.
def make_operon_fasta(operon_file, result_folder = './recover/blast_query_files/', operon_considered = 'ALL', filter_file = 'NONE'):

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
        #print len(result.keys())
        #print result['NC_000913'].keys()
        #print result['NC_000913']['atpA']
        #for genome in genome_of_interest_list:
        #    organism = genome.split('/')[len(i.split('/')) -1].split('.')[0]
        #    print organism
        #    org_dict = return_genbank_dict(genome)
        #    organism_dict_for_recovery.update({organism: org_dict})
    else: # TO DO: put  code in here to parse the filter file and then enforce that an NC number be in this file before evaluation
        pass
        
    #print len(organism_dict_for_recovery)

    exact_match = []
    missing_list = []

    for operon_line in [i.strip() for i in open(operon_file).readlines()]:
        operon = operon_line.split('\t')[0]
        if operon_considered == 'ALL':
            gene_list = operon_line.split('\t')[1:]
            #print "all operons will be used (but this is not working yet)."
        elif operon_considered == operon:
            gene_list = operon_line.split('\t')[1:]
            #print gene_list
            #print organism_dict_for_recovery.keys()
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
                    
            #print "Considering %s operon only" % operon
        else:
            pass
            #print operon, 'has been ignored'
    

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

# This function blasts every created DB against a query file.
# TO DO: Make this whole thing parallel
def blast(db_list, query_file, blast_result_folder, num_processors, eval_threshold):
    """Overview: This function performs a BLAST search of the supplied database and the list of query files dupplied.
       Then it saves the result in the gl_blast_folder.
       Returns: The result of the BLAST searches in the gl_blast_folder.
       Default: 
    """
    
    for db in db_list:
        out_file = "%s%s_prot.txt" % (blast_result_folder, db.split('/')[len(db.split('/')) - 1].split('.')[0])
        cmd = "blastall -p tblastn -a %i -i %s -d %s -e %s -o %s -m 9" % (num_processors, query_file, db, eval_threshold, out_file)
        os.system( cmd )
        

def perform_blast(db_folder, query_file, blast_result_folder, num_processors, eval_threshold):
    #print return_recursive_dir_files(db_folder)
    tmp = return_recursive_dir_files(db_folder)
    db_list = []
    for item in tmp:
        fname, file_ext = os.path.splitext(item)
        if file_ext == '.ffc':
            db_list.append(item)
        else:
            pass
        
    #print db_list
    
    blast(db_list, query_file, blast_result_folder, num_processors, eval_threshold)


# If we do not wish to filter any of the results, nc_filter_file should be set to the value 'NONE'
# I am only workign with the _prot right now, but we will have to deal with the rna thing too, later
def parallel_blast_parse(nc_filter_file, blast_folder, result_folder, operon_file, num_procs):
    if nc_filter_file == 'NONE':
        filtered_files = [(i, './recover/blast_parse/raw_info/')  for i in return_recursive_dir_files(blast_folder)]
    else:
        nc_filter = [i.strip() for i in open(nc_filter_file).readlines()]
        filtered_files =  [(i, './recover/blast_parse/raw_info/')  for i in return_recursive_dir_files(blast_folder) if i.split('/')[len(i.split('/'))-1].split('_p')[0] in nc_filter]
    #print len(filtered_files)
    #print filtered_files
    
    pool = Pool(processes = NUM_PROCESSORS)
    pool.map(initial_parse_blast, filtered_files)
    
    if nc_filter_file == 'NONE':
        filtered_files = [(i, result_folder) for i in return_recursive_dir_files(blast_folder)]
    else:
        nc_filter = [i.strip() for i in open(nc_filter_file).readlines()]
        filtered_files =  [(i, result_folder)  for i in return_recursive_dir_files(blast_folder) if i.split('/')[len(i.split('/'))-1].split('_p')[0] in nc_filter]
    
    pool = Pool(processes = NUM_PROCESSORS)
    pool.map(parse_blast, filtered_files)
    
    
# this is to be used with parallel_blast_parse ().  
# this is to make the first set of debugging files, sorry about it, but just have to deal with this while i
# am debugging things.
def initial_parse_blast(parse_tuple):
    fname, result_folder = parse_tuple
    hit_list = [i.strip() for i in open(fname).readlines() if i[0] != '#']

    outfile = './recover/blast_parse/raw_info/' + fname.split('/')[len(fname.split('/'))-1].split('_p')[0] + '.txt'

    handle = open(outfile, 'w')
    handle.write('\n'.join(hit_list))
    handle.close()
    
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
            
            
        tmp = item.split('\t')
        #NC_002516|Pseudomonas_aeruginosa_PAO1|PA5557|atpH|6252158|6252695|-1|Protein|NONE|60.7	ref|NC_000913|b3735|atpH|3917892|3918426|-1|Escherichia_coli_str._K-12_substr._MG1655	46.07	178	96	1	1	178	534	4	6e-51	 159
        ref = tmp[0].split('|')
        #print ref
        #print ref[0]
        #source_accession = ref[0].split('>')[1]
        source_accession = ref[0]
        source_common = ref[1]
        source_locus = ref[2]
        predicted_gene = ref[3]
        source_start = ref[4]
        source_start = ref[5]
        source_strand = ref[6]
        product_type = ref[7]
        synonyms = ref[8]
        synonym_list = ref[8].split(':')
        # fucking alert!!!! usign the wrong GC number here!!!!!
        source_gc = ref[9]
        accession, locus, gene, start, stop, strand, organism = tmp[1].split('|')[1:]
        percent_ident = tmp[2]
        e_val = tmp[10]
        bit_score = tmp[11]
        alignment_length = tmp[3]
        #operon_name = gene_to_operon_dict[predicted_gene]
        method = 'BLAST'
        if source_accession == accession:
            method = "exact"
        
        res = Homolog(accession, organism, locus, gene, predicted_gene, synonyms, e_val, percent_ident, bit_score, source_gc, start, stop, strand, product_type, alignment_length, method, source_accession, source_common, source_locus, source_start)    
        
        if locus not in homolog_locus_hit_dict.keys():
            homolog_locus_hit_dict.update({locus: []})
        homolog_locus_hit_dict[locus].append(res)
    print "homolog_locus_hit_dict", len(homolog_locus_hit_dict)

    
    result = []
    for locus in locus_hit_dict.keys():
        result.append(locus)
        for line in locus_hit_dict[locus]:
            result.append(line)
            
    
    
    outfile = result_folder + fname.split('/')[len(fname.split('/'))-1].split('_p')[0] + '.txt'
    print outfile
    outfile1 = result_folder + fname.split('/')[len(fname.split('/'))-1].split('_p')[0] + '_homolog.txt'
    print outfile1
    handle = open(outfile, 'w')
    handle.write('\n'.join(result))
    handle.close()
    
    
    print "got here"
    
    #hit_list = [i.strip() for i in open(fname).readlines() if i[0] != '#']
    
    '''for item in hit_list:
        locus = item.split('\t')[1].split('|')[2]
        if locus not in homolog_locus_hit_dict.keys():
            homolog_locus_hit_dict.update({locus: []})
        homolog_locus_hit_dict[locus].append(blast_hit_to_homolog(item))
        #if locus in homolog_locus_hit_dict.keys():
        #    homolog_locus_hit_dict[locus].append(blast_hit_to_homolog(item))
        #else:
        #    homolog_locus_hit_dict.update({locus: [blast_hit_to_homolog(item)]})'''

    
    homolog_result = []
    
    for locus in homolog_locus_hit_dict.keys():
        homolog_result.append(locus)
        for hlog in homolog_locus_hit_dict[locus]:
            homolog_result.append(hlog.ret_str())
            
    
    
    outfile1 = result_folder + fname.split('/')[len(fname.split('/'))-1].split('_p')[0] + '_homolog.txt'
    print outfile1
    handle1 = open(outfile1, 'w')
    handle1.write('\n'.join(homolog_result))
    handle1.close()
    
    
    

# Parse all blast results and per operon make a file containing all homologs, ordered alphabetically by NC number, saved in the result_folder.
# Returns a dictionary {operon:{organism:[homologs]}} 
# I might have to alter this... will have to see about all that..... ugh more thigns to do
def parse_blast_result(operon_file, blast_folder, homolog_folder, operon_folder, hgt_dict):
    file_list = os.listdir(blast_folder)
    gene_to_operon_dict = {}
    result = {}
    operon_list = [i.strip().split('\t') for i in open(operon_file).readlines()]
    #print operon_list
    for operon in operon_list:
        for gene in operon[1:]:
            gene_to_operon_dict.update({gene:operon[0]})
        result.update({operon[0]:{}})
    for fname in file_list:
       hit_list = [i.strip() for i in open(blast_folder + fname).readlines() if i[0] != '#']
       homolog_info = []
       for hit in hit_list:
           tmp = hit.split('\t')
           #NC_011833|Buchnera_aphidicola_str._5A_(Acyrthosiphon_pisum)|BUAP5A_006|atpA|4529|6068|1|Protein|NONE|30.6|exact	ref|NC_000907.1|HI0481|atpA|502547|504089|-1
           ref = tmp[0].split('|')
           source_accession = ref[0].split('>')[1]
           source_common = ref[1]
           source_locus = ref[2]
           source_start = ref[4]
           predicted_gene = ref[3]
           product_type = ref[7]
           synonyms = ref[8]
           synonym_list = ref[8].split(':')
           accession, organism, locus, gene, start, stop, strand, gc = tmp[1].split('|')
           percent_ident = tmp[2]
           e_val = tmp[10]
           bit_score = tmp[11]
           alignment_length = tmp[3]
           operon_name = gene_to_operon_dict[predicted_gene]
           method = 'blast'
           if source_accession == accession:
               method = "exact"
           #print accession
           #print hgt_dict.keys()
           '''if accession in hgt_dict.keys():
               if start in hgt_dict[accession].keys():
                   #hgt_status = {'likelyhood':'predicted_hgt','method':'alien-hunter','eval_score':float(c_val),'eval_thresh':float(c_thresh)}
                   hgt_status = hgt_dict[accession][start] #"likelyhood:predicted_hgt,method:alien-hunter,eval_score:%s,eval_thresh:%s" % (c_val, c_thresh)
                   homolog_info.append([accession, organism, locus, gene, predicted_gene, synonyms, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type, hgt_status])
                   if accession in result[operon_name].keys():
                        lst = result[operon_name].pop(accession)
                        lst.append(Homolog(accession, organism, locus, gene, predicted_gene, synonym_list, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type, hgt_status))
                        result[operon_name].update({accession:lst})
                   else: 
                       result[operon_name].update({accession:[Homolog(accession, organism, locus, gene, predicted_gene, synonym_list, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type, hgt_status)]})
               else:
                   homolog_info.append([accession, organism, locus, gene, predicted_gene, synonyms, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type])
                   if accession in result[operon_name].keys():
                       lst = result[operon_name].pop(accession)
                       lst.append(Homolog(accession, organism, locus, gene, predicted_gene, synonym_list, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type))
                       result[operon_name].update({accession:lst})
                   else: 
                       result[operon_name].update({accession:[Homolog(accession, organism, locus, gene, predicted_gene, synonym_list, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type)]})
           else:
               if accession in result[operon_name].keys():
                   lst = result[operon_name].pop(accession)
                   lst.append(Homolog(accession, organism, locus, gene, predicted_gene, synonym_list, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type))
                   result[operon_name].update({accession:lst})
               else: 
                   result[operon_name].update({accession:[Homolog(accession, organism, locus, gene, predicted_gene, synonym_list, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type)]})
                   
       handle = open(homolog_folder + fname, 'w')
       handle.write(gl_new_line.join(['\t'.join(i) for i in homolog_info]) + gl_new_line)
    for operon in result.keys():
        handle = open("%s%s_hits.txt" % (operon_folder, operon), 'w')
        for prokaryote in sorted(result[operon].keys()):
            handle.write(gl_new_line.join([i.ReturnHomologStr() for i in sorted(result[operon][prokaryote]) ]) + gl_new_line) '''
    return result


def blast_hit_to_homolog(hit):
    tmp = hit.split('\t')
    #NC_002516|Pseudomonas_aeruginosa_PAO1|PA5557|atpH|6252158|6252695|-1|Protein|NONE|60.7	ref|NC_000913|b3735|atpH|3917892|3918426|-1|Escherichia_coli_str._K-12_substr._MG1655	46.07	178	96	1	1	178	534	4	6e-51	 159
    ref = tmp[0].split('|')
    source_accession = ref[0].split('>')[1]
    source_common = ref[1]
    source_locus = ref[2]
    source_start = ref[4]
    source_start = ref[5]
    predicted_gene = ref[3]
    product_type = ref[7]
    synonyms = ref[8]
    synonym_list = ref[8].split(':')
    try:
        accession, organism, locus, gene, start, stop, strand, gc = tmp[1].split('|')
    except:
        print hit
    percent_ident = tmp[2]
    e_val = tmp[10]
    bit_score = tmp[11]
    alignment_length = tmp[3]
    operon_name = gene_to_operon_dict[predicted_gene]
    method = 'BLAST'
    if source_accession == accession:
        method = "exact"
        
    res = Homolog(accession, organism, locus, gene, predicted_gene, synonyms, e_val, percent_ident, bit_score, gc, start, stop, strand, product_type, alignment_length, method, source_accession, source_common, source_locus, source_start)
    return res


def main():

    # This file contains the phylogenetic distances from E. coli, as determined from rpoD. (unless it is rpoB). 
    # Either way, it contains two values. NC number and distance from E. coli, tab delineated.
    distance_file = './recover/Dist_test.txt'

    start = time.time()
    
    # this function should be called once before you start running the rest. not every time. ugh.
    #startup_crap(False)
    #startup_crap(True)
    
    #make_bast_db('./recover/AllGenbankPath.txt', NUM_PROCESSORS)
    
    closest_relatives_dict = return_closest_relatives(distance_file)

    #make_operon_fasta('./recover/operon_name_and_genes.txt', './recover/blast_query_files/', 'atpIBEFHAGDC')
     
    eval_threshold = .001
    
    #blast(db_list, query_files, blast_folder, NUM_PROCESSORS, eval_threshold)
    
    #perform_blast('./recover/db_folder/', './recover/blast_query_files/exact_matches.fa', './recover/blast_result/', NUM_PROCESSORS, eval_threshold)
    #parallel_blast('./recover/db_folder/', './recover/blast_query_files/exact_matches.fa', './recover/blast_result/', NUM_PROCESSORS, eval_threshold)
    
    
    #parse_blast_result()
    parallel_blast_parse('phylo_order.txt', './recover/blast_result/', './recover/blast_parse/intermediate/', './recover/operon_name_and_genes.txt', NUM_PROCESSORS)
    #parallel_blast_parse('NONE', './recover/blast_result/', './recover/blast_parse/intermediate/', './recover/operon_name_and_genes.txt', NUM_PROCESSORS)

    end = time.time()
    print end - start




if __name__ == '__main__':
    main()

