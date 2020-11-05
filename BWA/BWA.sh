cp $1 "bwa_genome1.fasta"
cp $2 "bwa_genome2.fasta"

printf "Starting...\n"

START_TIME=$(date +%s)

bwa index -a "is" "bwa_genome1.fasta"
printf "First file indexed\n\n"

bwa index -a "is" "bwa_genome2.fasta"
printf "Second file indexed\n\n"

bwa mem -o "bwa_output.sam" "bwa_genome1.fasta" "bwa_genome2.fasta"
# bwa bwasw -N 1 -f "bwa_output.sam" "bwa_genome1.fasta" "bwa_genome2.fasta"
printf "Aligned (.sam)\n\n"

# samtools view -Sb "bwa_output.sam" > "bwa_output.bam"
# printf "Converted to .bam\n\n"

# samtools sort -o "bwa_output_sorted.bam" "bwa_output.bam"
# printf "Sorted .bam\n\n"

# samtools index "bwa_output_sorted.bam"
# printf "Indexed .bam\n\n"

printf "Process time: $(( $(date +%s) - START_TIME ))s\n"

# sam2pairwise < bwa_output.sam > bwa_output_pairwise.txt

printf "Cleaning up...\n\n"

declare -a redundant_files=("bwa_genome1.fasta.amb"
			                "bwa_genome1.fasta.ann"
			                "bwa_genome1.fasta.bwt"
			                "bwa_genome1.fasta.pac"
			                "bwa_genome1.fasta.sa"

			                "bwa_genome2.fasta.amb"
			                "bwa_genome2.fasta.ann"
			                "bwa_genome2.fasta.bwt"
			                "bwa_genome2.fasta.pac"
			                "bwa_genome2.fasta.sa"

			                "bwa_genome1.fasta"
			                "bwa_genome2.fasta"
			                )

for file in "${redundant_files[@]}"
do
	if [[ -e "$file" ]]
	then
		rm -f "$file"
	fi
done
