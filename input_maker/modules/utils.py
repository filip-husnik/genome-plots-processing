# -*- coding: utf-8 -*-

import argparse
import csv
from Bio.Blast.Applications import NcbitblastxCommandline
from Bio import SeqFeature, SeqIO
import re
from time import localtime, strftime
from pathlib import Path
import sys



def current_time():
    """Returns the current time. When this function was executed."""
    return str(strftime("%Y-%m-%d %H:%M:%S", localtime()))


def get_args():

    parser = argparse.ArgumentParser(description="Create the necessary files to Visualize genome alignments",epilog="--------------------------")
    parser.add_argument("--order-file", help="Optional input file specifying order of the genomes in the figure. Genome record names should be listed on separate lines in order from top to bottom (as it will appear in the figure).", default=False)
    parser.add_argument("--evalue", default="1e-5", help="evalue for blast searches. Default is %(default)s.")

    args = parser.parse_args()
    
    return args


    
    
## Input: an integer or symbol
## Return type: string
## Use:Given an integer, output +,-,0 depending on it's value. Or N/A if not an int.
##     (For converting Bipython strand output to a symbol).
def strand(int1):
    feature_strand =  "None"
    
    if type(int1) == int:
        if int1 > 0:
            feature_strand = "+"
        elif int1 < 0:
            feature_strand = "-"
        else:
            feature_strand = "0"
    else:
        feature_strand = "None"
    
    return feature_strand

## Input: A biopython SeqFeature
## Return Type: Boolean
## Use: reports True if Feature type is not "gene" (not incl. pseudogenes) or "source"
def type_filter(seqFeature1):
    
    if type(seqFeature1) == SeqFeature.SeqFeature:
        
        if seqFeature1.type != "gene" and seqFeature1.type != "source" :
            return True
        elif "pseudo" in seqFeature1.qualifiers:
            return True
        else:
            return False
    else:
        return False

## Input: a Biopython SeqFeature
## Return type: String
## Use: return a string of the SeqFeature type 
##      Differs from SeqFeature.type because it includes 'pseudogene' as a type
def feature_type(seqFeature1):
    
    feature_type = "N/A"
    
    if type(seqFeature1) == SeqFeature.SeqFeature:
        
        if seqFeature1.type != "gene":
            feature_type = seqFeature1.type
        else:
            if "pseudo" in seqFeature1.qualifiers:
                feature_type = "pseudogene"
            else:
                feature_type = seqFeature1.type

    
    return feature_type


## Input: a Biopython SeqFeature
## Return type: String
## Use: return a string of the SeqFeature gene name if one exists.
##      Use locus_tag otherwise.
def get_gene(seqFeature1):
    
    gene_name = "N/A"
    
    if type(seqFeature1) == SeqFeature.SeqFeature:
        
        if "gene" in seqFeature1.qualifiers:
            gene_name = str(seqFeature1.qualifiers["gene"][0])
        else:
            if "locus_tag" in seqFeature1.qualifiers:
                gene_name = str(seqFeature1.qualifiers["locus_tag"][0])
            
    return gene_name


## Input: a Biopython SeqFeature
## Return type: Integer
## Use: returns an Integer of SeqFeature start position shifted by shift
##      Use locus_tag otherwise.
def start_position_shift(seqFeature1, shift):
    
    feature_start = -100
    
    if type(seqFeature1) == SeqFeature.SeqFeature:
        feature_start = seqFeature1.location.start + shift
            
    return feature_start

def get_dir_records(directory_name):

    if directory_name[-1] != "/":
        directory_name += "/"
    record_names=[]        # record names from genbank files
    genome_sizes = {}     #genome sizes in base pairs from genbank files
    records = []
    entries = Path(directory_name)
    for entry in entries.iterdir():
        gb_file = directory_name + entry.name
        for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):  #In case of multiple records
            record_names.append(gb_record.name.strip())
            genome_sizes[gb_record.name.strip()] = len(gb_record.seq)
            records.append(gb_record)
    
    return records, record_names, genome_sizes
    

def genbank_to_fasta(genbank_record, out_fasta):
    
    x = re.split("\/", out_fasta)
    x = x[-1]
    print('%s\tCreating %s' % (current_time(), x)),
    sys.stdout.flush()
    
    with open(out_fasta, "w") as f:
        SeqIO.write(genbank_record, f, "fasta")
    
    


def create_gene_info_csv(gb_record, out_csv):
    
    record_names=[]        # record names from genbank files
    genome_sizes = {}     #genome sizes in base pairs from genbank files
  
    recordName = gb_record.name.strip()
    record_names.append(recordName)
    genome_sizes[recordName] = len(gb_record.seq)
    
    x = re.split("\/", out_csv)
    x = x[-1]
    print('%s\tCreating gene info file: %s' % (current_time(), x)),
    sys.stdout.flush()
    
    with open(out_csv, mode="w", newline="") as gene_file:  ##newline="" to remove blanks b/w lines
      
        gene_writer = csv.writer(gene_file, delimiter=',')
        for feature in gb_record.features:
            
            if type_filter(feature):    #if not gene or source
                ##get strand symbol
                feature_strand = strand(feature.strand)
                ## get feature type. 
                featuretype = feature_type(feature)
                ## get gene name
                feature_gene = get_gene(feature)
                ## get feature start location, shifted by 1
                feature_start = start_position_shift(feature, 1)
                        
                gene_writer.writerow([feature_start, feature.location.end, featuretype, feature_strand, feature_gene])
    return record_names, genome_sizes

def get_genome_order(record_names, order_file):
    genome_order = []
    if order_file:
        with open(order_file, "r") as f:
            order_list = [line.strip() for line in f]
            for record in order_list:
                if record in record_names:
                    genome_order.append(record)
    else:
        genome_order = record_names
    
    return genome_order

def create_connections_csv(input_tsv,out_csv):
    
    x = re.split("\/", out_csv)
    x = x[-1]
    print('%s\tCreating connections file: %s' % (current_time(), x)),
    sys.stdout.flush()
    
    with open(out_csv, mode="w", newline="") as connections_file:
        connections_writer = csv.writer(connections_file, delimiter=',')
        with open(input_tsv, mode="r") as tsvfile:
            reader =  csv.reader(tsvfile, delimiter='\t')
            sign1 = "0"
            sign2 = "0"
            for row in reader:
                if int(row[7]) > 0:
                    sign1 = "+"
                elif int(row[7]) < 0:
                    sign1 = "-"
                if int(row[8]) > 0:
                    sign2 = "+"
                elif int(row[8]) < 0:
                    sign2 = "-"
                connections_writer.writerow([row[2], row[3], sign1, row[4], row[5], sign2])


    



def run_tblastx(args, query_fasta, subject_fasta, out_tsv):
    """Run BLASTX with query FASTA file and subject FASTA file."""
    
    x = re.split("\/", query_fasta)
    x = x[-1]
    y = re.split("\/", subject_fasta)
    y = y[-1]

    print('%s\ttBlastX executed with %s evalue cutoff. Query: %s. Subject: %s.' % (current_time(), args.evalue, x, y)),
    sys.stdout.flush()

    tblastx_cline = NcbitblastxCommandline(query=query_fasta,
                                         subject=subject_fasta,
                                         evalue=args.evalue,
                                         outfmt="\'6 qseqid sseqid qstart qend "
                                                "sstart send evalue qframe sframe\'",
                                         out=out_tsv)
    tblastx_cline()

    
    
    
    
    
