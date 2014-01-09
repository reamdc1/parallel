#!/usr/bin/python

#from multiprocessing import Pool
#import time
import os
import simplejson as json
import argparse

# Copyright(C) 2013 David Ream
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment


#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result

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

# This function will download the HGT files in islandviewer from Fiona Brinkman. This is just a quick thing to update the files if I need this
# to run later.

def download_islandviewer_hgt(dest_folder):
    
    # This file specifically does not tell me what method was used to determine the data, so it is not very informative for the purpose we want.
    #download_remote_file('http://www.pathogenomics.sfu.ca/islandviewer/download_all.php?tool=islandviewer&type=csv', './hgt/all_gis_islandviewer.csv')

    download_remote_file('http://www.pathogenomics.sfu.ca/islandviewer/download_all.php?tool=islandpick&type=csv', "%sall_gis_islandpick.csv" % dest_folder)
    
    download_remote_file('http://www.pathogenomics.sfu.ca/islandviewer/download_all.php?tool=islandpath_dimob&type=csv', "%sall_gis_islandpath_dimob.csv" % dest_folder)
    
    download_remote_file('http://www.pathogenomics.sfu.ca/islandviewer/download_all.php?tool=sigi_hmm&type=csv', "%sall_gis_sigi_hmm.csv" % dest_folder)
    

# return 
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

# Currently I am not doing anything with this piece of code, since I think that it should be used i the main code fork.
# This piece of code is supposed to assign HGT to orfs in an organism.  I suspect that it requires a bit more work.
def assign_hgt(operon_folder, hgt_dict):
    flist = sorted(returnRecursiveDirFiles(operon_folder))
    summary_char = ['$$', '@@', '##']
    result = []
    for operon in flist:
        print operon
        lines = [i.strip() for i in open(operon).readlines() if i[:2] not in summary_char]
        for line in lines:
            if len(line) <= 2: # there are sometimes oddities in the files that I wish to ignore
                pass
            elif line[:2] != '++':
                NC, org, locus, annotation, predicted, syn, e_val, a, b, c, start, stop, strand, seq_type, hgt_stuff = line.split('\t')
                start = int(start)
                stop = int(stop)
                if NC in hgt_dict.keys():
                    hgt_list = hgt_dict[NC]
                    for item in hgt_list:
                        hgt_start, hgt_stop, hgt_size, hgt_program = item
                        if start >= hgt_start and start <= hgt_stop:
                            print NC, org, predicted, hgt_program, 'hgt suspected!' 
                        else:
                            pass
                        
                else:
                    pass
            else:
                pass

def main():

    parser = argparse.ArgumentParser(description='Download precomputed HGT files or parse existing ones into a python dictionary that is stored in JSON.')

    parser.add_argument("-f", "--folder", dest="folder", metavar="FOLDER", default='./',
                help="Folder where results will be downloaded or read from, depending on the mode selected. Default is the folder where the script is run from.")
                
    parser.add_argument("-o", "--outfile", dest="outfile", default='hgt.json', metavar="FILE",
                help="Location where the program will store the json formatted result of the program. If omitted the program will use the folder plus 'HGT.json'.")
                
    parser.add_argument("-d", "--download", dest="download", action="store_true", default=False,
                help="Add this option if you wish to download a new set of HGT files from the IslandPath website.")

    parsed_args = parser.parse_args()
    folder = parsed_args.folder
    outfile = folder + parsed_args.outfile
    download = parsed_args.download
    
    if download:
        download_islandviewer_hgt(folder)
    
    
    flist = returnRecursiveDirFiles(folder)
    hgt_dict = parse_hgt_files(flist)
    
    handle = open(outfile, 'w')
    json.dump(hgt_dict, handle)
    handle.close()

if __name__ == '__main__':
    main()
