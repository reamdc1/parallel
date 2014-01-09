#!/usr/bin/env python
import sys
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import os

# take the genbank file specified by genbank path, and save the customized result file in the db_directory folder

# I am modifying this program so that i can 

def convert_genbank(genbank_path, db_directory):

    record_list = []
    seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
    nc_number = seq_record.name
    
    # genbank information about the sequence found in the file. 
    accession = seq_record.id
    comment = seq_record.annotations['comment']
    data_file_division = seq_record.annotations['data_file_division'] # this is not very useful, since it is historical and does not reflect current taxonomy
    date = seq_record.annotations['date'] # this should be the date that the file was updated last. uesful if i wish to autoupdate for more current versions of a genome
    gi = seq_record.annotations['gi'] # i do not plan on using this, but it may be handly to have, so i will include anyway
    key_words = seq_record.annotations['keywords']
    organism = seq_record.annotations['organism'].replace(' ', '_')
    sequence_version = seq_record.annotations['sequence_version'] # this is an integer
    taxonomy = seq_record.annotations['taxonomy'] # this will be a list
    
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

    print 'organism ', organism
    
    
    
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
