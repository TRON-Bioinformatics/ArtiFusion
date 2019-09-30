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

## Usage

### Input Data
The tool requires the following files as input data:
 - Input table with fusion partners and size ratio thresholds (Header/Columns: Gene symbol 1;Ratio;Gene symbol 2)
 - Gene Model in BED format (Example can also be found in test_folder/test.bed)
 - HGNC gene symbol mapping table (Example can also be found in test_folder/test_gene_symbols.csv)
 - Reference Genome as fasta (Can be downloaded from http://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)

### Prepare the references
The input reference files can be generated as follows:
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

### Output

ArtiFuse produces a summary file with the generated ArtiFusions as well as a modified FASTA file.

 - `summary.csv` - Summary of the generated ArtiFusions including gene symbols, transcript names, breakpoint positions and the replacement sequence
 - `simulated.fa` - Modified reference genome assembly, mostly used for fusion detection tool index generation
 - `<chrom>.fa` - For each chromosome a file with the according modified sequence is generated (useful for MapSplice index generation)

Table 1 - Columns description for file `summary.csv`

| **Column** | **Description** |
|:-----------|:----------------|
| **Gene\_Symbol\_1** | Gene symbol of the 5' end fusion partner |
| **Gene\_Symbol\_2** | Gene symbol of the 3' end fusion partner |
| **Transcript\_1** | Transcript ID of the 5' end fusion partner |
| **Transcript\_2** | Transcript ID of the 3' end fusion partner |
| **BP1** | Chromosomal position of the 5' end of fusion junction; 1-based coordinate |
| **BP2** | Chromosomal position of the 3' end of fusion junction; 1-based coordinate |
| **Exp\_Ratio** | Ratio from Input Table |
| **Obs\_Ratio** | Actual Ratio after swapping |
| **Swap\_Sequence** | Replacement Sequence from Gene A being inserted in Gene B |