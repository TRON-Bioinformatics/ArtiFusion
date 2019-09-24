#!/usr/bin/env python


import sys

from operator import itemgetter
from argparse import ArgumentParser


class GTF2BED(object):
    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile

    def run(self):
        
        trans_dict = {}

        with open(self.infile) as inf:
            for line in inf:
                if line.startswith("#"):
                    continue

                elements = line.rstrip().split("\t")
                chrom = elements[0]
                type = elements[2]
                start = int(elements[3])
                end = int(elements[4])
                strand = elements[6]

                info = elements[8]

                info_dict = {}
                info_eles = info.rstrip(";").split("; ")
                for ele in info_eles:
                    key, val = ele.split(" ", 1)
                    info_dict[key] = val.replace("\"","")



                if type == "transcript":
                    trans_id = info_dict["transcript_id"]
                    trans_dict[trans_id] = {"chrom": chrom, "start": start-1, "end": end, "strand": strand, "exons": []}

                elif type == "exon":
                    trans_id = info_dict["transcript_id"]
                    trans_dict[trans_id]["exons"].append((start-trans_dict[trans_id]["start"]-1, end-start+1))


        outf = open(self.outfile, "w")
        

        for trans_id in sorted(trans_dict, key=lambda x: (trans_dict[x]["chrom"], trans_dict[x]["start"])):
            outf.write("\t".join([
                trans_dict[trans_id]["chrom"], 
                str(trans_dict[trans_id]["start"]), 
                str(trans_dict[trans_id]["end"]), 
                trans_id, 
                "0", 
                trans_dict[trans_id]["strand"], 
                str(trans_dict[trans_id]["start"]), 
                str(trans_dict[trans_id]["end"]), 
                "0", 
                str(len(trans_dict[trans_id]["exons"])), 
                ",".join(str(y) for x, y in sorted(trans_dict[trans_id]["exons"])) + ",",
                ",".join(str(x) for x, y in sorted(trans_dict[trans_id]["exons"])) + ",",
            ]) + "\n"
            )

        outf.close()

def main():
    parser = ArgumentParser(description='Converts GTF to BED format')
    parser.add_argument('-i', '--input', dest='input_file', help='Specify input GTF file', required=True)
    parser.add_argument('-o', '--output', dest='output_file', help='Specify output file', required=True)

    args = parser.parse_args()

    gb = GTF2BED(args.input_file, args.output_file)
    gb.run()


if __name__ == "__main__":
    main()
