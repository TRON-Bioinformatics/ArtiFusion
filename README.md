# ArtiFusion

ArtiFusion is a tool to simulate artificial fusion events by modifying a given reference genome.
The tool copies parts of the exonic sequence of gene A within the reference genome FASTA sequence into the downstream region of gene B and replaces the copied regions of gene A with Ns.
The breakpoints are defined by using a size ratio between gene A and gene B and are always placed on exon-exon junctions.
Intronic and intergenic regions remain unchanged.

The approach can be used to benchmark fusion detection tools with realistic biological data.
In contrast to simulating NGS reads (ART package, https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm), we do not lose the biological relevance of sequencing data.


## Installation

```
git clone https://github.com/TRON-Bioinformatics/ArtiFusion.git
```

###  Dependencies

 - Python 2.7.15
 - Python packages:
   - [Biopython](https://biopython.org/) (>= 1.72)

Data:
 - Input table with fusion partners and size ratio thresholds (Header/Columns: Gene symbol 1;Ratio;Gene symbol 2)
 - Gene Model in BED format (Example can also be found in test_folder/test.bed)
 - HGNC gene symbol mapping table (Example can also be found in test_folder/test_gene_symbols.csv)
 - Reference Genome as fasta (Can be downloaded from http://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)

## Usage

### Prepare the references
```
wget http://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.97.gtf.gz

python misc/gtf2bed.py -i Homo_sapiens.GRCh38.97.gtf -o <transcripts_bed>
python misc/gtf2symbol.py -i Homo_sapiens.GRCh38.97.gtf -o <transcript_to_genesymbol_tsv>

```
### Execute the tool
```
python artifuse.py \
  -i test_folder/test_input_table.csv \
  -g <primary_assembly_fasta> \
  -b <transcripts_bed> \
  -f <transcript_to_genesymbol_tsv> \
  -o <working_dir>
```
	      