#!/usr/bin/env python

import sys

from argparse import ArgumentParser

class GTF2Symbol(object):
    def __init__(self, gtf, out):
        self.gtf = gtf
        self.outfile = out

    def run(self):
        trans_to_gene = {}
        with open(self.gtf) as inf:
            for line in inf:
                elements = line.rstrip().split("\t")

                if len(elements) < 3:
                    continue

                if elements[2] == "transcript":
                    trans_id = ""
                    gene_symbol = ""
                    ids = elements[8].split(";")
                    for i in range(len(ids)):
                        if "transcript_id" in ids[i]:
                            trans_id=ids[i].rsplit(" ",1)[1]
                        elif "gene_name" in ids[i]:
                            gene_symbol=ids[i].rsplit(" ",1)[1]

                    trans_to_gene[trans_id.strip("\"")] = gene_symbol.strip("\"")

        outf = open(self.outfile, "w")

        for key in trans_to_gene.iterkeys():
            outf.write(key + "\t" + trans_to_gene[key] + "\n")

        outf.close()

def main():
    parser = ArgumentParser(description='Converts GTF to Gene Symbol mapping file')
    parser.add_argument('-i', '--input', dest='input_file', help='Specify input GTF file', required=True)
    parser.add_argument('-o', '--output', dest='output_file', help='Specify output file', required=True)

    args = parser.parse_args()

    gs = GTF2Symbol(args.input_file, args.output_file)
    gs.run()

                
if __name__ == '__main__':
    main()
