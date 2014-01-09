#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import simplejson as json
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Copyright(C) 2013 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment


NUM_PROCESSORS = os.sysconf("SC_NPROCESSORS_CONF")

# convert a file of format : operon name then a list of the full names of the genes within that operon
# into a dictionary that can be easily accessed or filtered later.
def parse_operon_file(fname):
    result = {}
    
    for line in [i.strip().split('\t') for i in open(fname).readlines()]:
        result.update({line[0]: line[1:]})
        
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
                    header = '|'.join([accession, common_name, locus, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
                    result.update({new_gene: [header, ''.join(seq)]})
                    unk_cnt +=1
                else:
                    header = '|'.join([accession, common_name, locus, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
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
def make_operon_fasta(gene_list, genbank_list, num_processors, folder, refrence):

    pool = Pool(processes = num_processors)
    organism_dict_for_recovery = dict(pool.map(parallel_operon_fasta, genbank_list))
    
    protein_match = []
    rna_match = []
    pseudogene_match = []
    missing_list = []
    
    refrence_prot = []
    refrence_rna = []

    # This list should be updated should other types of RNAs be annotated later as parts of (siRNA comes to mind).
    RNA_codes = ['rRNA', 'tRNA', 'ncRNA']
    print RNA_codes
    for org in organism_dict_for_recovery.keys():
        print "org", org
        for gene in gene_list:
            if gene in organism_dict_for_recovery[org].keys():
                print organism_dict_for_recovery[org][gene][0]
                if organism_dict_for_recovery[org][gene][0].split('|')[7] == 'Protein':
                    outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                    id=organism_dict_for_recovery[org][gene][0], description = '')
                    protein_match.append(outseq)
                    #if org == refrence:
                    #    outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                    #    id=organism_dict_for_recovery[org][gene][0], description = '')
                    #    refrence_prot.append(outseq)
                elif organism_dict_for_recovery[org][gene][0].split('|')[7] in RNA_codes:
                    print "RNA"
                    outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                    id=organism_dict_for_recovery[org][gene][0], description = '')
                    rna_match.append(outseq)
                    #if org == refrence:
                    #    print org, refrence
                    #    outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                    #    id=organism_dict_for_recovery[org][gene][0], description = '')
                    #    refrence_rna.append(outseq)
                    #else:
                    #    print "Fail", org, refrence
                elif organism_dict_for_recovery[org][gene][0].split('|')[7] == 'Pseudo_Gene':
                    outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                    id=organism_dict_for_recovery[org][gene][0], description = '')
                    pseudogene_match.append(outseq)
                    #if org == refrence:
                    #    outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                    #    id=organism_dict_for_recovery[org][gene][0], description = '')
                    #    refrence_prot.append(outseq)
                else:
                    print organism_dict_for_recovery[org][gene][0]
            else: # The gene is missing, we will look for it at a later stage in the program
                #print organism_dict_for_recovery[org][gene][0]
                item = '\t'.join([org, gene])
                missing_list.append(item)
            

        
        handle = open(folder + 'protein_matches.fa', 'w')       
        SeqIO.write(protein_match, handle,"fasta")
        handle.close()
        
        handle = open(folder + 'rna_matches.fa', 'w')       
        SeqIO.write(rna_match, handle,"fasta")
        handle.close()
        
        handle = open(folder + 'pseudogene_matches.fa', 'w')       
        SeqIO.write(pseudogene_match, handle,"fasta")
        handle.close()
        
        handle = open(folder + 'missing_operon_genes.txt', 'w')
        handle.write('\n'.join(missing_list))
        handle.close()
        
        handle = open(folder + 'refrence_prot.fa', 'w')
        SeqIO.write(refrence_prot, handle,"fasta")
        handle.close()
        
        handle = open(folder + 'refrence_rna.fa', 'w')
        SeqIO.write(refrence_rna, handle,"fasta")
        handle.close()
        
        
    
        
# At a later time should we deem it necessary, this is where we should place all the code that will test all of command line params
# and make sure that they are sensible.
def option_check(parsed_args):
    pass



def categorize_operons(ref_org, genbank_list, operon_dict):
    ref_path = [i for i in genbank_list if i.split('/')[-1].split('.')[0] == ref_org][0]
    
    prot_operons = []
    mixed_operons = []
    ref_dict = return_genbank_dict(ref_path)
    
    #print ref_dict
    
    for operon in sorted(operon_dict.keys()):
        print operon
        for gene in operon_dict[operon]:
            try:
                gene_type = ref_dict[gene][0].split('|')[7]
                if gene_type == 'Protein' or gene_type == 'Pseudogene':
                    print "Protein", gene, gene_type
                else:
                    print "RNA", gene, gene_type
                
                #print gene, gene_type, ref_dict[gene]
            except:
                print "Missing gene %s found in refrence organism." % gene

def download_remote_file(url, outfile, block = 4096): # low memory implementation, if you have big files.
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


def regulon_db(outfile):

    url = 'http://regulondb.ccg.unam.mx/menu/download/datasets/files/OperonSet.txt'
    download_remote_file(url, "./OperonSet.txt")
    
    # this bit of code parses the regulon DB file, keep in mind the columns:
    # (1) Operon name
    # (2) First gene-position left
    # (3) Last gene-position right
    # (4) DNA strand where the operon is coded
    # (5) Number of genes contained in the operon
    # (6) Name or Blattner number of the gene(s) contained in the operon
    # (7) Evidence that support the existence of the operon's TUs
    
    tmp_list = []
    for line in [i.strip() for i in open("./OperonSet.txt").readlines()]:
        if len(line) < 2:
            pass
        elif line[0] == '#':
            pass
        else:
            tmp_list.append(line.strip())
    #print tmp_list
    
    result = []
    for line in tmp_list:
        #print line.split('\t')
        try:
            name, gene_pos_left, gene_pos_right, strand, length, gene_names, evidence = line.split('\t')
        except:
            name, gene_pos_left, gene_pos_right, strand, length, gene_names = line.split('\t')
        if int(length) >=5:
            print name, gene_names
            result.append([name] + gene_names.split(','))
    print result
    
    handle = open(outfile, 'w')
    handle.write('\n'.join(['\t'.join(i) for i in result]))
    handle.close()
            
    


def main():
    
    start = time.time()
    
    parser = argparse.ArgumentParser(description='Make several fasta files containing all annotated operon genes from every organism that is under consideration for BLAST searching. Files created are protein_matches.fa,rna_matches.fa, pseudogene_matches.fa, and missing_operon_genes.txt which cannot be changed currently.')

    parser.add_argument("-i", "--infile", dest="infile", metavar="FILE", default='./operon_name_and_genes.txt',
                help="Input operon file. It should be parsed into tab delienated format with operon name first, then each gene contained.")
                
    parser.add_argument("-o", "--operon", dest="operon_name", default='', metavar="STRING",
                help="Location where the program will store the fasta. If omitted the program will use the name 'operon_multiple_refrence.fa'.")
                
    parser.add_argument("-g", "--genbank", dest="genbank_file_list", metavar="FILE", default = './genbank_pathway_lists/filtered_list_genbank_paths.txt',
                help="File that contains the full pathway to every genbank file you wish to consider.")
                
    parser.add_argument("-f", "--folder", dest="folder", metavar="FOLDER", default = './',
                help="Folder where all output will be stored. By default results are stored in the current directory.")
                
    parser.add_argument("-r", "--refrence", dest="refrence", metavar="STRING", default = 'NC_000913',
                help="Accession number of the refrence organism. This information is used to determine the product type of each gene (RNA/Protein), a necessary piece of information to classify the operons that are under investigation.")
                
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")

    parser.add_argument("-d", "--download", dest="download", action="store_true", default=False,
                help="Add this option if you wish to download the regulonDB operon file and parse it, instead of using something older.")
    
    parsed_args = parser.parse_args()
    infile = parsed_args.infile
    folder = parsed_args.folder
    ref_org = parsed_args.refrence
    download_regulondb = parsed_args.download
    
    
    if download_regulondb:
        regulon_db(infile)
    
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = parsed_args.num_proc
    
    genbank_list = [i.strip() for i in open(parsed_args.genbank_file_list).readlines()]
    
    operon_dict = parse_operon_file(infile)
    
    # Code block that handles the targeting of a single operon.
    gene_list = []
    if parsed_args.operon_name == '':
        for operon in operon_dict.keys():
            gene_list = gene_list + operon_dict[operon]
    else:
        try:
            gene_list = gene_list + operon_dict[parsed_args.operon_name]
        except:
            print "The operon: %s is not found." % parsed_args.operon_name
            sys.exit()
    
    print "genbank_list ", genbank_list        
    make_operon_fasta(gene_list, genbank_list, num_proc, folder, ref_org)
    
    
    categorize_operons(ref_org, genbank_list, operon_dict)
    
    # easy way to run this, using all the defaults that make sense
    # ./make_operons_fasta.py -f ./fasta_generation/
    print time.time() - start
    
if __name__ == '__main__':
    main()
