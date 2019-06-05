# ArtiFusion

ArtiFusion is a tool to simulate artificial fusions in the reference genome


## Install

```
git clone https://github.com/TRON-Bioinformatics/ArtiFusion.git
```

###  Dependencies

Python packages:

 - [Biopython](https://biopython.org/) (>= 1.72)

Data:

 - HGNC gene symbol mapping table
 - Reference genome as fasta

## Usage

```
python artifuse.py \
  -i test_folder/<input_table> \
  -g test_folder/test.fa \
  -b test_folder/test.bed \
  -f test_folder/test_gene_symbols.csv \
  -o test_out
```
	      