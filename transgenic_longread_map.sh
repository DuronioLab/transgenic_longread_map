#!/bin/bash

#SBATCH -p general
#SBATCH --ntasks=2
#SBATCH --time=5:00:00
#SBATCH --mem=16g

#################################################################################################################################
#                                                                                                                               #
#   Run the script with:                                                                                                        #
#   sbatch --time=5:00:00 --mem=16g --ntasks=2 --wrap="sh ./scripts/alignment_script.sh <ref.fasta> <query.fasta> [min_length]" #
#                                                                                                                               #
#################################################################################################################################

# Check if at least two FASTAs were provided
if [[ -z "$2" ]]; then
    echo "Usage: $0 <reference_fasta_file> <query_fasta_file> [length]"
    exit 1
fi

if [[ $1 != *.fa && $1 != *.fasta ]] || [[ $2 != *.fa && $2 != *.fasta ]]; then
    echo "Error: The files must end with .fa or .fasta"
    exit 1
fi

# Grab the reference and query FASTAs
ref_fasta=$1
query_fasta_file=$2

# The length is the third argument, or 5000 if not provided
length=${3:-5000}

# Print the filenames and length
echo "You have provided:"
echo "Found FASTA reference file: "${ref_fasta}
echo "Query Filename: "${query_fasta_file}
echo "The size of the cutoff is: "${length}


#Extract all .gz fastq files
gunzip *.fastq.gz

## Concatenate all FASTQ files together
# Find all files generated from the GridION ("FAX"..) in the current directory with the ".fastq" extension
files=$(find . -maxdepth 1 -name "FAX*.fastq")

# Concatenate the found files
cat $files > concat.fastq

# Remove the file extension from the file name
ref_basename=$(basename "$ref_fasta")
ref_basename="$(ref_basename%.*}"

echo "The base that will be used for naming is " ${ref_basename}
echo

#Rename each FASTQ read and set up query.fasta
printf "\n Renaming FASTQ reads and generating FASTA file\n"
module purge && module load seqkit

seqkit replace --quiet -p .+ -r "seq_{nr}" concat.fastq > renamed_reads.fastq

seqkit seq --quiet --min-len $length renamed_reads.fastq -o out.fastq
seqkit fq2fa --quiet out.fastq -o out.fasta

printf "\nAfter filtering there are this many reads:\n"
echo $(cat out.fastq|wc -l)/4|bc

split --lines=40000 -d --additional-suffix=".fasta" out.fasta out.

#Search for reference query
module purge && module load blat
printf "\nPerforming BLAT and restarting sequences\n"

# Initialize an empty file to store all results
touch all_results.psl

# Loop through each fasta file
for fasta_file in out.*.fasta
do
    # Run the blat command on the current fasta file
    blat "$fasta_file" query.fasta -oneOff=3 -noHead output.psl

    # Append the results to the all_results.psl file
    cat output.psl >> all_results.psl
done

awk 'NR > 6 {print $14}' all_results.psl > temp_names.txt
sort temp_names.txt | uniq > blat_names.txt

#For each read, restart based on blat
module load seqkit

readarray -t File < blat_names.txt
touch filtered_reads.fasta
for f in "${File[@]}"
do
  echo ${f} > fname.txt
  seqkit grep --quiet -n -f fname.txt out.fasta >> filtered_reads.fasta
done

#Convert back to FASTQ with fake quality scores
module purge && module load seqtk
seqtk seq -F '#' filtered_reads.fasta > ${ref_basename}_filtered.fastq

module load minimap2

minimap2 --secondary=no --sam-hit-only -ax map-ont ${ref_fasta} ./${ref_basename}_filtered.fastq > ./${ref_basename}_mapped.sam

module purge && module load samtools

# Sort BAM file
samtools sort ./${ref_basename}_mapped.sam -o ./${ref_basename}_mapped.bam
samtools view -bq 1 ./${ref_basename}_mapped.bam > ./${ref_basename}_unique.bam


#Generate the consensus sequence
module purge && module load medaka
printf "\nPerforming Medaka consensus search\n"
out_folder="${ref_basename}_consensus"

printf "\nMedaka will use the reference FASTA: \n"${ref_fasta}"\n"

medaka_consensus -i ./filtered_reads.fasta -o ${ref_basename} -d ${ref_fasta} -t 1
