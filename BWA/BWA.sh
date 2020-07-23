cp $1 "bwa_genome1.fasta"
cp $2 "bwa_genome2.fasta"

printf "Starting...\n"

START_TIME=$(date +%s)

bwa index "bwa_genome1.fasta"
printf "First file indexed\n\n"

bwa index "bwa_genome2.fasta"
printf "Second file indexed\n\n"

bwa mem -o "bwa_output.sam" "bwa_genome1.fasta" "bwa_genome2.fasta"
printf "Aligned (.sam)\n\n"

samtools view -Sb "bwa_output.sam" > "bwa_output.bam"
printf "Converted to .bam\n\n"

printf "Process time: $(( $(date +%s) - START_TIME ))s\n"

sam2pairwise < bwa_output.sam > bwa_output_pairwise.txt