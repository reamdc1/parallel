#!/usr/bin/python

import time
import os
import sys
import argparse
from multiprocessing import Pool
from Bio import SeqIO
from Bio.SeqUtils import GC

# Copyright(C) 2013 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment

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

def regulon_db(outfile, number_of_genes, url, download, experimantal):
    if download:
        download_remote_file(url, "./regulonDB/OperonSet.txt")
    
    # this bit of code parses the regulon DB file, keep in mind the columns:
    # (1) Operon name
    # (2) First gene-position left
    # (3) Last gene-position right
    # (4) DNA strand where the operon is coded
    # (5) Number of genes contained in the operon
    # (6) Name or Blattner number of the gene(s) contained in the operon
    # (7) Evidence that support the existence of the operon's TUs
    
    tmp_list = []
    for line in [i.strip() for i in open("./regulonDB/OperonSet.txt").readlines()]:
        if len(line) < 2:
            pass
        elif line[0] == '#':
            pass
        else:
            tmp_list.append(line.strip())

    result = []
    for line in tmp_list:
        try:
            name, gene_pos_left, gene_pos_right, strand, length, gene_names, evidence = line.split('\t')
            #print len(evidence.split('xperiment'))
        except:
            name, gene_pos_left, gene_pos_right, strand, length, gene_names = line.split('\t')
        if int(length) >=number_of_genes and experimantal:
            #print name, gene_names
            if len(evidence.split('xperiment')) > 1:
                result.append([name] + gene_names.split(','))
        elif int(length) >=number_of_genes and not experimantal:
            result.append([name] + gene_names.split(','))

    handle = open(outfile, 'w')
    handle.write('\n'.join(['\t'.join(i) for i in result]))
    handle.close()

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
            #print synonyms
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
                    #################################################
                    # New code to overcome some issues with         #
                    # RegulonDB. we will store synonym data as well #
                    # which improves operon recovery slightly.      #
                    #################################################
                    if synonyms != 'NONE':
                        #print synonyms
                        for syn in synonyms.split(':'):
                            #print syn
                            header = '|'.join([accession, common_name, locus, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
                            result.update({syn: [header, ''.join(seq)]})

    #print 'The number of non-protein regions in %s is: %i.' % (common_name, cnt)
    return result

# This function will allow me to do the main work of make_operon_fasta, but allow parallel
# processing. I will have to make a parallel array, which will take some time to learn.  i 
# will keep this stub for later to implement. 
def parallel_genome_dict(genome):
    organism = genome.split('/')[-1].split('.')[0]
    organism_dict_for_recovery = {}
    org_dict = return_genbank_dict(genome)
    organism_dict_for_recovery.update({organism: org_dict})
    return (organism, org_dict)
  
def main():
    
    start = time.time()

    parser = argparse.ArgumentParser(description="Download and parse a regulonDB operon file, then determine information about the component genes by investigating the refrence organisms.")
                
    parser.add_argument("-o", "--outfile", dest="outfile", metavar="FILE", default='./regulonDB/operon_names_and_genes.txt',
                help="Output of the initial RegulonDB file.")
    
    parser.add_argument("-i", "--infile", dest="infile", default='./refrence_organism_paths.txt', metavar="FILE",
                help="A file that contains the path to every refrence organism to be used (The file I use has both sequenced varients of the E.coli K-12 genome).")
                
    parser.add_argument("-f", "--folder", dest="folder", default='./regulonDB/', metavar="FOLDER",
                help="A folder that will contain all of the output files that are renerated with this program.")
    
    parser.add_argument("-n", "--num_genes", dest="min_num_genes", metavar="INT", default = 5, type=int,
                help="Minum number of genes contained within an operon to consider it part of the dataset.")
    
    parser.add_argument("-u", "--url", dest="url", default='http://regulondb.ccg.unam.mx/menu/download/datasets/files/OperonSet.txt', metavar="URL",
                help="A file that contains the BLAST query for every gene of interest in the dataset.")
                
    parser.add_argument("-d", "--download", dest="download", action="store_false", default=True,
                help="Add this option if you wish to download the regulonDB operon, otherwise the program will assume that you have already done this step.")

    parser.add_argument("-e", "--experimantal", dest="experimental_only", action="store_false", default=True,
                help="Add this option if you wish to download the regulonDB operon, otherwise the program will assume that you have already done this step.")
    
    parsed_args = parser.parse_args()
    outfile = parsed_args.outfile
    infile = parsed_args.infile
    min_number_of_genes = parsed_args.min_num_genes
    url = parsed_args.url
    download = parsed_args.download
    experimental = parsed_args.experimental_only
    folder = parsed_args.folder
    
    genbank_list = [i.strip() for i in open(infile).readlines()]
    
    pool = Pool(processes = 2)
    organism_dict_for_recovery = dict(pool.map(parallel_genome_dict, genbank_list))
    
    #print organism_dict_for_recovery.keys()

    regulon_db(outfile, min_number_of_genes, url, download, experimental)
    
    protein_only_list = []
    mixed_list = [] # this is a list that contains operons that contain both protein coding and RNA coding genes

    gene_dict = {}
    
    # open the file that contains the pathways to each of the refrence genomes (since this is regulonDB, both varients of E. coli apply
    for operon_line in [i.strip().split('\t') for i in open(outfile).readlines()]:
        operon_name = operon_line[0]
        gene_dict.update({operon_name:{}})
        #print operon_name
        for gene in operon_line[1:]:
            gene_dict[operon_name].update({gene:[]})
            for ref_org in sorted(organism_dict_for_recovery.keys()):
                try:
                    gene_product =  organism_dict_for_recovery[ref_org][gene][0].split('|')[7]
                    gene_dict[operon_name][gene].append(gene_product)
                except:
                    pass

    # there are two classes of potential products, protein and rna.  there is a list of terms that are acceptable for each
    # These two lists keep track of this, so we know what are effectively congruent terms, as we only care if they prot/rna are signaled.
    # when we determine the typ of each operon, which is RNA, protein, mixed.
    rna_list = ['tRNA', 'rRNA', 'ncRNA', ]
    protein_list = ['Protein', 'Pseudo_Gene']
    
    prot_result = []
    rna_result = []
    mixed_result = []
    
    # this can be removed later, after I later change downstream operon file parsing
    old_result = []

    for operon in sorted(gene_dict.keys()):
        operon_error = False # If there is missing/inconsistent information in the genes then ignore the operon
        operon_type = ''
        gene_list = []
        old_gene_list = []
        for gene in gene_dict[operon]:
            gene_type = ''
            if len(gene_dict[operon][gene]) == 0:
                #print operon, gene, "Missing annotation information"
                operon_error = True
            elif len(gene_dict[operon][gene]) == 1:
                if gene_dict[operon][gene][0] in protein_list:
                    gene_type = 'p'
                elif gene_dict[operon][gene][0] in rna_list:
                    gene_type = 'r'
                else:
                    pass
            else:
                tmp_name = gene_dict[operon][gene][0]
                for next_type in gene_dict[operon][gene][1:]:
                    if next_type in protein_list and tmp_name in protein_list:
                        gene_type = 'p'
                    elif next_type in rna_list and tmp_name in rna_list:
                        gene_type = 'r'
                    else:
                        #print operon, gene, "Annotation information disagrees"
                        operon_error = True


            gene_list.append("%s:%s" % (gene, gene_type))
            old_gene_list.append(gene)
            if operon_error:
                pass
            elif operon_type == '': # We have not evaluated the opern type (protein, RNA, mixed) yet, so set it to the type of the first gene
                operon_type = gene_type
            elif operon_type != gene_type:
                operon_type = 'm'
            elif operon_type == gene_type:
                pass
            else:
                print "the function broke"
        if operon_error:
            print "rejected operon ", operon
            operon_error = False
        elif operon_type == 'p':
            prot_result.append('\t'.join([operon] + gene_list))
            # Remove next line after operon list parsing has been updated
            old_result.append('\t'.join([operon] + old_gene_list))
        elif operon_type == 'r':
            rna_result.append('\t'.join([operon] + gene_list))
        else:
            mixed_result.append('\t'.join([operon] + gene_list))
    #print "prot_result", prot_result
    #print "rna_result", rna_result
    #print "mixed_result", mixed_result
    
    
    
    
    #############################################################################################################################################
    # currently I am dumping 4 files out into the out folder. They are not selectable in terms of names, sorry. I may fix this, i may not.      #
    # the new operon format, where i have the gene names and the type of product that they code for is more useful. also we have the ability    #
    # to select prot only, rna only, or both.  To include every operon, we still would have to cat the files, which i have not done, but i      #
    # think might be a good idea.  I will look into that once we are planning on conisdering all operons. This is a simple operation. I am not  #
    # going to do it now, because after 3+ days sosrting out all the bugs in this program, files, format, etc... i cannot stand the thought     #
    # of working further on this program.                                                                                                       #
    #############################################################################################################################################

    handle_prot = open(folder + 'operon_name_and_genes_prot_only.txt', 'w')
    handle_rna = open(folder + 'operon_name_and_genes_rna_only.txt', 'w')
    handle_mixed_type = open(folder + 'operon_name_and_genes_mixed_type.txt', 'w')
    
    handle_prot.write('\n'.join(prot_result))
    handle_prot.close()
    
    handle_rna.write('\n'.join(rna_result))
    handle_rna.close()
    
    handle_mixed_type.write('\n'.join(mixed_result))
    handle_mixed_type.close()
    
    handle = open(folder+ "operon_name_and_genes.txt", 'w')
    handle.write('\n'.join(old_result))
    handle.close()
    
    
    
    
    print time.time() - start
   
if __name__ == '__main__':
    main()