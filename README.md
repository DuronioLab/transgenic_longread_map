# transgenic_longread_map

[![DOI](https://zenodo.org/badge/791369050.svg)](https://zenodo.org/doi/10.5281/zenodo.11060915)

Script for use in UNC Longleaf that will take whole genome long read data (FASTQ), a reference (FASTA), and a set of unique sequences to anchor the long reads in (FASTA) and will generate an alignment (SAM/BAM) and a consensus sequence.

# Running

Run with
```
sbatch --time=5:00:00 --mem=16g --ntasks=2 --wrap="sh ./scripts/alignment_script.sh <ref.fasta> <query.fasta> [min_length]"
```

# Requirements
- sbatch (though you can run tools separately or without sbatch)
- blat
- seqkit
- seqtkk
- minimap2
- medaka
- samtools

# Required files

## Ref.fasta
A well assembled reference of any size. A reference of just the transgenic locus with the expected insertions is recommended. Errors in this reference can adversely affect the accuracy.

## query.fasta
A FASTA file with 50-70 bases of sequence denoting unique sequence that will be used to anchor any reads that have ambiguous sequence. For example: [loxP]-[BAC]-[4xDWT] 50 bases of the "loxP element are sufficient to anchor any reads to this region if another [BAC]-[4xDWT] is present. This also prevents endogenous histone genes from aligning to a transgenic locus. Tested with up to 8 query sequences, but more may be desirable.

# Optional arguments

## min_length
The minimum length of reads that will be kept prior to alignment. Default is set to 5000 bp.
