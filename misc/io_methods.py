#!/usr/bin/env python

import os
import pwd
import grp
import stat

from Bio import SeqIO

def create_folder(path):
    '''This function creates a specified folder, if not already existing and grants the respective permission.'''
    if not os.path.exists(path):
        print("Creating folder", path)
        os.makedirs(path)

def grant_permissions(path, mode):
    '''This function grants specific permissions for a given folder or file.'''
    for root, dirs, files in os.walk(path, topdown=False):
        for dir in dirs:
            os.chmod(dir, mode)
        for file in files:
            os.chmod(file, mode)

def get_fastq_files(paths):
    '''This function returns a list of fastq files for the given list of paths.'''
    fastqs = []
    for path in paths:
        if os.path.isdir(path):
            files = os.listdir(path)
            for file in files:
                file_path = os.path.join(path,file)
                if os.path.isfile(file_path) and file.endswith((".fastq.gz","fastq")):
                    fastqs.append(file_path)
                elif os.path.isdir(file_path):
                    fastqs_tmp = IOMethods.get_fastq_files([file_path])
                    fastqs.extend(fastqs_tmp)
        elif os.path.isfile(path) and path.endswith((".fastq.gz","fastq")):
            fastqs.append(path)
    return fastqs

def get_gene_symbol_map(infile):
    '''This function loads the reference gene symbol mapping.'''
    print("Loading gene symbols into memory... (input={})".format(infile))
    gene_map = {}
    with open(infile) as inf:
        for line in inf:
            elements = line.rstrip().split("\t")
            id = elements[0]
            gene_symbol = elements[1]
            gene_map.setdefault(gene_symbol, []).append(id)
    return gene_map

def get_trans_map(infile):
    '''This function loads the reference transcripts database into memory.'''
    print("Loading transcriptome into memory... (input={})".format(infile))
    trans_map = {}
    with open(infile) as inf:
        for line in inf:
            elements = line.rstrip().split("\t")
            chrom = elements[0]
            start = int(elements[1])
            end = int(elements[2])
            id = elements[3]
            strand = elements[5]
            block_sizes = [int(x) for x in elements[10].rstrip(",").split(",")]
            block_starts = [int(x) for x in elements[11].rstrip(",").split(",")]
            abs_block_starts = []
            for i, exon in enumerate(block_starts):
                abs_block_starts.append(start + block_starts[i])
            trans_map[id] = [chrom, start, end, strand, block_starts, abs_block_starts, block_sizes]
    return trans_map


def get_fasta_map(infile):
    '''This function loads the genome assembly fasta into memory.'''
    print("Loading genome into memory... (input={})".format(infile))
    chr_map = {}
    for seq_record in SeqIO.parse(infile, 'fasta'):
        chr_map[seq_record.id] = str(seq_record.seq)
#    chr_map = {}
#    chrom = ""
#    seq_str = ""
#    seq = []
#    with open(infile) as inf:
#        for line in inf:
#            if line.startswith(">"):
#                if chrom != "":
#                    chr_map[chrom] = seq
#                    seq_str = "".join(seq)
#                    chr_map[chrom] = seq_str
#                    print("Chromosome {} loaded into memory. ({}bp)".format(chrom, len(seq_str)))

#                chrom = line.rstrip()[1:]
#                seq = []
#                seq = ""
#            else:
#                seq.append(line.rstrip())
#                seq += line.rstrip()
#        chr_map[chrom] = seq
#        seq_str = "".join(seq)
#        chr_map[chrom] = seq_str
#        print("Chromosome {} loaded into memory. ({}bp)".format(chrom, len(seq_str)))
    return chr_map
