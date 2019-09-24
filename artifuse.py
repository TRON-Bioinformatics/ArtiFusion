#!/usr/bin/env python

import os
import sys
import string

from argparse import ArgumentParser

import misc.io_methods as IOMethods

class ArtiFusion(object):
    def __init__(self, input_table, working_dir, reference_genome, bed_file, gene_symbols_file):
        self.input_table = input_table
        self.working_dir = working_dir
        self.swap_pos = []
        self.compl_dict = string.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx")

        self.gene_map = IOMethods.get_gene_symbol_map(gene_symbols_file)
        self.trans_map = IOMethods.get_trans_map(bed_file)
        self.chr_map = IOMethods.get_fasta_map(reference_genome)

    def output_fasta(self):
        '''This function outputs the modified reference genome in FASTA format.'''
        outf = open(os.path.join(self.working_dir, "simulated.fa"), "w")
        for key in sorted(self.chr_map):
            outf_split = open(os.path.join(self.working_dir, "{}.fa".format(key)), "w")
            outf.write(">{}\n".format(key))
            outf_split.write(">{}\n".format(key))
            chrom_seq = self.chr_map[key]
            for i in range(0, len(chrom_seq), 60):
                outf.write(chrom_seq[i:i + 60] + "\n")
                outf_split.write(chrom_seq[i:i + 60] + "\n")
            print("Modified chromosome {} written to file. ({}bp)".format(key, len(chrom_seq)))
            outf_split.close()
        outf.close()

    def get_longest_transcript(self, gene_symbol):
        '''This function fetches the longest transcript of a gene from the BED genes table.'''
        max_size = 0
        sel_id = ""
        ids = self.gene_map[gene_symbol]

        for id in ids:
            trans_size = sum(self.trans_map[id][6])
            if trans_size > max_size:
                max_size = trans_size
                sel_id = id

        return sel_id

    def rev_comp(self,seq):
        '''This function generates the reverse complement of a transcriptomic sequence.'''
        return seq.translate(self.compl_dict)[::-1]


    def determine_bp_1(self, abs_block_starts, block_sizes, chrom, ratio, size, strand):
        '''This function returns the breakpoint on Gene 1 depending on strandedness.
        Insertion sequence is in italic. Arrow marks strandedness. BP1=Breakpoint 1
        If the strand is positive this returns the end of the leftmost exon:
        |>>>>>>>>|<BP1>-----|>>>>>>>>>>|------|>>>>>>>>>>>|
        If the strand is negative this returns the start of rightmost exon:
        //////////---------//////////--<BP1>||||||||||||<--
        '''
        swap_seq = ""
        swap_size = 0
        bp = ""
        swaps = []
        if strand == "+":
            for i in range(len(abs_block_starts)-1, 0, -1):
                exon_start = abs_block_starts[i]
                exon_end = exon_start + block_sizes[i]
                exon_seq = self.chr_map[chrom][exon_start:exon_end]

                print("Exon {}: ({},{})".format(i, exon_start, exon_end))
                print("{} : {}".format(len(swap_seq + exon_seq), (ratio*size)))
                swap_seq = exon_seq + swap_seq
                swap_size += block_sizes[i]
                swaps.append((chrom,exon_start,exon_end))

                bp = "{}:{}:{}".format(chrom, abs_block_starts[i-1] + block_sizes[i-1], strand)
                if len(swap_seq) > ratio * size:

                    break
        elif strand == "-":
            for i in range(0, len(abs_block_starts)-1, 1):
                exon_start = abs_block_starts[i]
                exon_end = exon_start + block_sizes[i]
                exon_seq = self.chr_map[chrom][exon_start:exon_end]

                print("Exon {}: ({},{})".format(i,exon_start,exon_end))
                print("{} : {}".format(len(swap_seq + exon_seq), (ratio*size)))

                swap_seq += exon_seq
                swap_size += block_sizes[i]
                swaps.append((chrom,exon_start,exon_end))


                bp = "{}:{}:{}".format(chrom, abs_block_starts[i+1], strand)
                if len(swap_seq) > ratio * size:
                    break
            swap_seq = self.rev_comp(swap_seq)
        return (swap_seq, swap_size, bp, swaps)

    def determine_bp_2(self, abs_block_starts, block_sizes, chrom, swap_size, strand):
        '''This function returns the breakpoint on Gene 2 depending on strandedness.
        Sequence to be replaced is in italic. Arrow marks strandedness. BP2=Breakpoint 2
        If the strand is positive this returns the end of the leftmost exon:
        -->||||||||||->>-|||||||||||-->>-<BP2>//////////-->>--////////////
        If the strand is negative this returns the start of rightmost exon:
        //////////--<<--//////////<BP2>-<<--||||||||||||-<<--||||||||||<--
        '''
        swap_seq = ""
        swap_exon_start = 0
        bp = ""
        if strand == "+":
            for i in range(len(abs_block_starts) - 1, 0, -1):
                exon_start = abs_block_starts[i]
                exon_end = exon_start + block_sizes[i]
                exon_seq = self.chr_map[chrom][exon_start:exon_end]

                swap_seq += exon_seq
                if len(swap_seq) > swap_size:
                    bp = "{}:{}:{}".format(chrom, exon_start + 1, strand)
                    swap_exon_start = i
                    break
        elif strand == "-":
            for i in range(len(abs_block_starts)):
                exon_start = abs_block_starts[i]
                exon_end = exon_start + block_sizes[i]
                exon_seq = self.chr_map[chrom][exon_start:exon_end]

                swap_seq += exon_seq
                if len(swap_seq) > swap_size:
                    bp = "{}:{}:{}".format(chrom, exon_end + 1, strand)
                    swap_exon_start = i
                    break

        return (swap_exon_start, bp)


    def swap(self, abs_block_starts, block_sizes, chr_2, strand_2, swap_exon_start, swap_seq, swap_size):
        '''This function replaces the region within Gene 2 with the region from Gene 1.'''
        pos = 0
        swapped = False
        replacement_seq = ""
        if strand_2 == "+":
            abs_pos = 0
            for i in range(swap_exon_start, len(abs_block_starts)):
                    
                exon_start = abs_block_starts[i]
                exon_end = exon_start + block_sizes[i]
                replacement_seq += self.chr_map[chr_2][abs_pos:exon_start]
                abs_pos = exon_end
                if swapped:
                    replacement_seq += "N" * block_sizes[i]
                else:
                    if pos + block_sizes[i] <= swap_size:
                        swap_region = swap_seq[pos:pos+block_sizes[i]]
                        orig_region = self.chr_map[chr_2][exon_start:exon_end]
                        replacement_seq += swap_region
                    else:
                        cutoff = pos + block_sizes[i] - swap_size
                        swap_region = swap_seq[pos:swap_size]
                        orig_region = self.chr_map[chr_2][exon_start:exon_end-cutoff]
                        replacement_seq += swap_region + "N"*cutoff
                        swapped = True
                pos += block_sizes[i]
            replacement_seq += self.chr_map[chr_2][abs_pos:]
        elif strand_2 == "-":
            abs_pos = len(self.chr_map[chr_2])
            for i in range(swap_exon_start,-1,-1):
                exon_start = abs_block_starts[i]
                exon_end = exon_start + block_sizes[i]
                replacement_seq = self.chr_map[chr_2][exon_end:abs_pos] + replacement_seq
                abs_pos = exon_start
                if swapped:
                    replacement_seq = "N" * block_sizes[i] + replacement_seq
                else:
                    if pos + block_sizes[i] <= swap_size:
                        swap_region = swap_seq[pos:pos+block_sizes[i]]
                        replacement_seq = self.rev_comp(swap_region) + replacement_seq
                    else:
                        cutoff = pos + block_sizes[i] - swap_size
                        swap_region = swap_seq[pos:swap_size]
                        replacement_seq = "N" * cutoff + self.rev_comp(swap_region) + replacement_seq
                        swapped = True
                pos += block_sizes[i]

            replacement_seq = self.chr_map[chr_2][:abs_pos] + replacement_seq
        return replacement_seq

    def run(self):
        '''This function iterates over the input table and determines BP positions and generates a modified reference genome FASTA.'''
        outf_summary = open(os.path.join(self.working_dir, "summary.csv"), "w")
        outf_summary.write("Gene_Symbol_1;Transcript_1;BP1;Gene_Symbol_2;Transcript_2;BP2;Exp_Ratio;Obs_Ratio;Swap_Sequence\n")

        with open(self.input_table) as inf:
            next(inf)
            for line in inf:
                elements = line.rstrip().split(";")
                gene_symbol_1 = elements[0]
                ratio = float(elements[1])
                gene_symbol_2 = elements[2]
                swap_exon_start = 0
                

                sel_id_1 = ""
                sel_id_2 = ""
                try:                
                    sel_id_1 = self.get_longest_transcript(gene_symbol_1)
                except:
                    continue

                (chr_1, start, end, strand, block_starts, abs_block_starts, block_sizes) = self.trans_map[sel_id_1]
                size = sum(block_sizes)
                
                (swap_seq, swap_size, bp_1, swaps) = self.determine_bp_1(abs_block_starts, block_sizes, chr_1, ratio, size, strand)
                attempts = 0
                while swap_size == 0 and attempts < 40:
                    attempts += 1
                    ratio += 0.01
                    (swap_seq, swap_size, bp_1, swaps) = self.determine_bp_1(abs_block_starts, block_sizes, chr_1, ratio, size, strand)

                if swap_size == 0 or swap_size == size:
                    print("Swap size still zero after attempting to increase gene 1 ratio...")
                    outf_summary.write("{}_{}: Starting exon of gene 1 cannot be 1. Gene ratio (Current: {}).\n".format(gene_symbol_1, gene_symbol_2, ratio))
                    continue
                elif swap_size > 0 and swap_size != size and attempts > 0:
                    print("Ratio increase successful. New ratio: {}".format(ratio))
                self.swap_pos.extend(swaps)
                try:
                    sel_id_2 = self.get_longest_transcript(gene_symbol_2)
                except:
                    continue
                (chr_2, start, end, strand, block_starts, abs_block_starts, block_sizes) = self.trans_map[sel_id_2]
                size = sum(block_sizes)

                (swap_exon_start, bp_2) = self.determine_bp_2(abs_block_starts, block_sizes, chr_2, swap_size, strand)

                if swap_exon_start == 0:
                    print("Swap Exon Start zero... (Ratio: {})".format(ratio))
                    outf_summary.write("{}_{}: Starting exon of gene 2 cannot be 1. Modify gene ratio (Current: {}).\n".format(gene_symbol_1, gene_symbol_2, ratio))
                    continue

                self.chr_map[chr_2] = self.swap(abs_block_starts, block_sizes, chr_2, strand, swap_exon_start, swap_seq, swap_size)

                print("{};{};{};{};{};{};{};{};{}".format(gene_symbol_1,
                                                          sel_id_1,
                                                          bp_1,
                                                          gene_symbol_2,
                                                          sel_id_2,
                                                          bp_2,
                                                          ratio,
                                                          len(swap_seq) / size,
                                                          swap_seq))
                outf_summary.write("{};{};{};{};{};{};{};{};{}\n".format(gene_symbol_1,
                                                                         sel_id_1,
                                                                         bp_1,
                                                                         gene_symbol_2,
                                                                         sel_id_2,
                                                                         bp_2,
                                                                         ratio,
                                                                         len(swap_seq) / size,
                                                                         swap_seq))


        outf_summary.close()

        # Replace Exons left from BP with N's
        for (chrom, start, end) in self.swap_pos:
            self.chr_map[chrom] = self.chr_map[chrom][0:start] + "N" * (end - start) + self.chr_map[chrom][end:]

        IOMethods.create_folder(self.working_dir)
        self.output_fasta()


def main():
    parser = ArgumentParser(description="Simulates fusion genes based on wildtype genes")

    parser.add_argument("-i", "--input-table", dest="input_table", help="Specify input file with genes to simulate fusions.", required=True)
    parser.add_argument("-g", "--reference-genome", dest="reference_genome", help="Specify reference genome assembly (.fa) to modify.", required=True)
    parser.add_argument("-b", "--bed", dest="bed_file", help="Specify reference genes file.", required=True)
    parser.add_argument("-f", "--gene-symbols", dest="gene_symbols_file", help="Specify reference gene symbols file.", required=True)
    parser.add_argument("-o", "--working-dir", dest="working_dir", help="Specify path to output FASTA file containing modified genome.", required=True)
    

    args = parser.parse_args()

    af = ArtiFusion(args.input_table, args.working_dir, args.reference_genome, args.bed_file, args.gene_symbols_file)

    af.run()

if __name__ == '__main__':
    main()
