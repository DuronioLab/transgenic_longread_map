#!/bin/bash

#SBATCH -p general
#SBATCH --ntasks=2
#SBATCH --time=5:00:00
#SBATCH --mem=16g

#################################################################################################################################################################################################################
#                                                                                                                                                                                                               #
#   Run the script with:                                                                                                                                                                                        #
#   sbatch --time=5:00:00 --mem=16g --ntasks=2 --wrap="sh ./scripts/transgenic_longread_map.sh --ref <ref.fasta> --query <query.fasta> -- size [min_length] --fastqs <.FASTQ, .FASTQ.GZ, directories, or mix>"  #
#                                                                                                                                                                                                               #
#######################################################################################################################################################################1#########################################
# Initialize our own variables
ref_fasta=""
query_fasta_file=""
length=0
found_fastq=""
declare -a fastqs

# A function to display a usage message
usage() { echo "Usage: $0 [--ref REF] [--query QUERY] [--size SIZE] [--fastqs FASTQS]" 1>&2; exit 1; }


# Process the arguments
while (( "$#" )); do
  case "$1" in
    --ref)
      ref_fasta="$2"
      shift 2
      ;;
    --query)
      query_fasta_file="$2"
      shift 2
      ;;
    --size)
      length="$2"
      shift 2
      ;;
    --fastqs)
      while [[ "$2" != --* && -n "$2" ]]; do
        fastqs+=("$2")
        shift
      done
      shift
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done

# Process --fastqs
for i in "${fastqs[@]}"; do
    # process "$i"
    if [[ -d $i ]]; then
        # If directory, find all .fastq files recursively
        fastq_files=$(find $i -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) ! -name "concat.fastq")
        #echo "Fastq files in directory $i:"
        #echo "$fastq_files"
    elif [[ -f $i ]]; then
        # If file, check if it ends with .fastq and is not named concat.fastq
        if [[ ($i == *.fastq || $i == *.fastq.gz) && $i != "concat.fastq" ]]; then
            #echo "Fastq file: $i"
            fastq_files=$i
        fi
    else
        echo "$i is not a valid file or directory"
    fi
    found_fastq+=$fastq_files$'\n'
done

if ! [[ $length =~ ^[0-9]+$ ]] || [ $length -le 0 ]; then
    echo "Warning: --size should be an integer greater than 0. Defaulting to 5000"
    length=5000
fi

if [[ $ref_fasta != *.fa && $ref_fasta != *.fasta ]]; then
  echo "Error --ref should be a file that ends in .fa or .fasta"
  exit 1
fi

if [[ $query_fasta_file != *.fa && $query_fasta_file != *.fasta ]]; then
  echo "Error --query should be a file that ends in .fa or .fasta"
  exit 1
fi

# Print the filenames and length
echo "You have provided:"
echo -e "Found the FASTQ files: \n"${found_fastq}
echo "Found FASTA reference file: "${ref_fasta}
echo "Query Filename: "${query_fasta_file}
echo "The size of the cutoff is: "${length}

# Create an empty file to store the concatenated results
cat /dev/null > concat.fastq

# Loop over all arguments
for file in ${found_fastq[@]}
do
    # Check if the file is a .fastq.gz file
    if [[ $file == *.fastq.gz ]]
    then
        # If it is, gunzip the file and append to the output file
        gunzip -c "$file" >> concat.fastq
    elif [[ $file == *.fastq ]]
    then
        # If it's a .fastq file, just append it to the output file
        cat "$file" >> concat.fastq
    else
        echo "Warning: $file is not a .fastq or .fastq.gz file and will be skipped."
    fi
done

module load nanopack
NanoPlot --fastq_rich concat.fastq --N50 -o ./NanoPlot

# Remove the file extension from the file name
ref_basename=$(basename "$ref_fasta")
ref_basename="${ref_basename%.*}"

echo "The base that will be used for naming is " ${ref_basename}
echo

total_reads=$(($(wc -l < concat.fastq) / 4))

if (( total_reads > 100000 ))
then
  proportion=$(echo "100000 / $total_reads" | bc -l)
  proportion=$(printf "%.2f" $proportion)
else
  proportion=1.0
fi

echo "Keeping ${proportion} of reads"

module load seqkit

# Randomly select 1000 reads with seqkit and process them
seqkit sample -p $proportion concat.fastq > concat_50k.fastq

awk '{
    if (NR % 4 == 1 || NR % 4 == 3) {
        print $0
    } else {
        seq_length=length($0)
        start=(length / 2) - 50
        print substr($0, start, 100)
    }
}' concat_50k.fastq > concat_50k_100bp.fastq


/proj/mckaylab/genomeFiles/fastq_screen_v0.11.1/fastq_screen --force --aligner bowtie2 -conf ./fastq_screen.conf concat_50k_100bp.fastq --outdir ./FQscreen/

rm concat_50k.fastq
rm concat_50k_100bp.fastq

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
minimap2 --secondary=no --sam-hit-only -ax map-ont dm6_mt_wMel.fasta ./concat.fastq > ./${ref_basename}_dm6mtwMel.sam

module purge && module load samtools

# Sort BAM file
samtools sort ./${ref_basename}_mapped.sam -o ./${ref_basename}_mapped.bam
samtools view -bq 1 ./${ref_basename}_mapped.bam > ./${ref_basename}_unique.bam

samtools sort ./${ref_basename}_dm6mtwMel.sam -o ./${ref_basename}_dm6mtwMel.bam
samtools index ./${ref_basename}_dm6mtwMel.bam

samtools idxstats ./${ref_basename}_dm6mtwMel.bam | cut -f 1,3 | grep -v "\*" > output.txt

echo "set terminal pdf size 8,6" > plot.gp
echo "set output '${ref_basename}_dm6mtwMel.pdf'" >> plot.gp
echo "set style data histograms" >> plot.gp
echo "set style fill solid 1.0 border -1" >> plot.gp
echo "set xtics rotate by -45" >> plot.gp
echo "plot 'output.txt' using 2:xtic(1) title 'Number of reads aligned to each chromosome'" >> plot.gp
gnuplot plot.gp
rm output.txt plot.gp

#Generate the consensus sequence
module purge && module load medaka
printf "\nPerforming Medaka consensus search\n"
out_folder="${ref_basename}_consensus"

printf "\nMedaka will use the reference FASTA: \n"${ref_fasta}"\n"

medaka_consensus -i ./filtered_reads.fasta -o ${ref_basename} -d ${ref_fasta} -t 1
