# -*- coding: utf-8 -*-

"""
Requires: Biopython installed
"""

import os
import csv
import glob
import re
import sys
from modules import utils

## assumes that every record has a seq name, as well as genome size in bps


args = utils.get_args()

script_path = os.path.dirname(os.path.realpath(__file__)) + "/"

genbank_folder = 'genbank_files/'
data_folder = 'data/'
gbfile_directory = script_path + genbank_folder #directory with genbank files
data_directory = script_path + data_folder            #directory to store data files


gb_records,record_names,genome_sizes = utils.get_dir_records(gbfile_directory) #get all records from a dir with genbank files

for gb_record in gb_records:
    out_path = data_directory + gb_record.name.strip()
    utils.genbank_to_fasta(gb_record, out_path + ".fasta")  #create fasta from gb files
    utils.create_gene_info_csv(gb_record, out_path + "_gene_info_processing.csv") #create gene info csv files

        

genome_order = utils.get_genome_order(record_names, args.order_file)


print('%s\tCreating data_files.csv' % (utils.current_time())),
sys.stdout.flush()

with open(data_directory + "data_files.csv", mode="w", newline="") as datafile:
    data_writer = csv.writer(datafile, delimiter=',')
    data_writer.writerow(["genome", "genes", "connections"])
    
    for i in range(len(genome_order) - 1):
        
        g1 = genome_order[i]
        g2 = genome_order[i+1]
        
        query_fasta = data_directory + g1 + ".fasta"
        subject_fasta = data_directory + g2 + ".fasta"
        out_tsv = data_directory + g1 + "_vs_" + g2 + ".tsv"
        out_csv = data_directory + g1 + "_vs_" + g2 + ".csv"
        
        utils.run_tblastx(args, query_fasta, subject_fasta, out_tsv) 
        utils.create_connections_csv(out_tsv,out_csv)
        
        data_writer.writerow([genome_sizes[g1], data_folder + g1 + "_gene_info_processing" + ".csv", data_folder + g1 + "_vs_" + g2 + ".csv" ])
    data_writer.writerow([genome_sizes[genome_order[-1]], data_folder + genome_order[-1] + "_gene_info_processing" + ".csv" ])


for f in glob.glob(data_directory + "*.fasta"):
    
    x = re.split("\/", f)
    x = x[-1]
    print('%s\tRemoving %s' % (utils.current_time(), x)),
    sys.stdout.flush()
    os.remove(f)

for f in glob.glob(data_directory + "*.tsv"):
    
    x = re.split("\/", f)
    x = x[-1]
    print('%s\tRemoving %s' % (utils.current_time(), x)),
    sys.stdout.flush()
    os.remove(f)

        
    
##TODO: Colour schemes, single or double strand?
##TODO: Create data directory?
    

         
