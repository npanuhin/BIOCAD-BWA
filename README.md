<h1 align="center">BWA</h1>

## Result

- Small:
	- deletion
		- [**SAM**](./BWA/small/deletion/bwa_output.sam)
		- [**BAM**](./BWA/small/deletion/bwa_output.bam)
		- [**pairwise**](./BWA/small/deletion/bwa_output_pairwise.txt) (View with disabled *word wrap*)
	- insertion
		- [**SAM**](./BWA/small/insertion/bwa_output.sam)
		- [**BAM**](./BWA/small/insertion/bwa_output.bam)
		- [**pairwise**](./BWA/small/insertion/bwa_output_pairwise.txt) (View with disabled *word wrap*)
	- inversion
		- [**SAM**](./BWA/small/inversion/bwa_output.sam)
		- [**BAM**](./BWA/small/inversion/bwa_output.bam)
		- [**pairwise**](./BWA/small/inversion/bwa_output_pairwise.txt) (View with disabled *word wrap*)
	- translocation
		- [**SAM**](./BWA/small/translocation/bwa_output.sam)
		- [**BAM**](./BWA/small/translocation/bwa_output.bam)
		- [**pairwise**](./BWA/small/translocation/bwa_output_pairwise.txt) (View with disabled *word wrap*)
- Large:
	- [**SAM**](./BWA/large/bwa_output.sam)
	- [**BAM**](./BWA/large/bwa_output.bam)
	- [**pairwise**](./BWA/large/bwa_output_pairwise.txt) (View with disabled *word wrap*)


## Contents

- `large/large_genome1.fasta`: [source](https://www.ncbi.nlm.nih.gov/nuccore/CP003305.1)
- `large/large_genome2.fasta`: [source](https://www.ncbi.nlm.nih.gov/nuccore/CP000766.3)

## Requirements

- **BWA**: http://bio-bwa.sourceforge.net
- **Samtools**: https://www.htslib.org
- **sam2pairwise**: https://github.com/mlafave/sam2pairwise

`sudo apt install bwa samtools`


## Run

- Navigate to `BWA/small/{...}` or `BWA/large` folder
- Run `./run.sh`

###### Pairwise

`sam2pairwise < bwa_output.sam > bwa_output_pairwise.txt`