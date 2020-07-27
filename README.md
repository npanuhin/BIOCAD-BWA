<h1 align="center">BWA</h1>

> This repository is not intended to represent work, but rather to store and transmit information.

## Result

For `SAM` and `pairwise` files *word wrap* should be disabled

- Small:

	[**source fasta**](./samples/small/source.fasta)

	- deletion
		- [**deletion fasta**](./samples/small/deletion.fasta)
		- [**SAM**](./BWA/small/deletion/bwa_output.sam)
		- [**BAM**](./BWA/small/deletion/bwa_output.bam)
		- [**pairwise**](./BWA/small/deletion/bwa_output_pairwise.txt)

	- insertion
		- [**insertion fasta**](./samples/small/insertion.fasta)
		- [**SAM**](./BWA/small/insertion/bwa_output.sam)
		- [**BAM**](./BWA/small/insertion/bwa_output.bam)
		- [**pairwise**](./BWA/small/insertion/bwa_output_pairwise.txt)

	- inversion
		- [**inversion fasta**](./samples/small/inversion.fasta)
		- [**SAM**](./BWA/small/inversion/bwa_output.sam)
		- [**BAM**](./BWA/small/inversion/bwa_output.bam)
		- [**pairwise**](./BWA/small/inversion/bwa_output_pairwise.txt)

	- inversion2
		- [**inversion2 fasta**](./samples/small/inversion2.fasta)
		- [**SAM**](./BWA/small/inversion2/bwa_output.sam)
		- [**BAM**](./BWA/small/inversion2/bwa_output.bam)
		- [**pairwise**](./BWA/small/inversion2/bwa_output_pairwise.txt)

	- translocation
		- [**translocation fasta**](./samples/small/translocation.fasta)
		- [**SAM**](./BWA/small/translocation/bwa_output.sam)
		- [**BAM**](./BWA/small/translocation/bwa_output.bam)
		- [**pairwise**](./BWA/small/translocation/bwa_output_pairwise.txt)

- Large1:
	- [**genome1 fasta**](./samples/large1/large_genome1.fasta)
	- [**genome2 fasta**](./samples/large1/large_genome2.fasta)
	- [**SAM**](./BWA/large1/bwa_output.sam)
	- [**BAM**](./BWA/large1/bwa_output.bam)
	- [**pairwise**](./BWA/large1/bwa_output_pairwise.txt)

- Large2:
	- [**genome1 fasta**](./samples/large2/large_genome1.fasta)
	- [**genome2 fasta**](./samples/large2/large_genome2.fasta)
	- [**SAM**](./BWA/large2/bwa_output.sam)
	- [**BAM**](./BWA/large2/bwa_output.bam)
	- [**pairwise**](./BWA/large2/bwa_output_pairwise.txt)

- Large3:
	- [**genome1 fasta**](./samples/large3/large_genome1.fasta)
	- [**genome2 fasta**](./samples/large3/large_genome2.fasta)
	- [**SAM**](./BWA/large3/bwa_output.sam)
	- [**BAM**](./BWA/large3/bwa_output.bam)
	- [**pairwise**](./BWA/large3/bwa_output_pairwise.txt)


## Stages

1. **BWA** indexes two `fasta` sequences
2. **BWA** aligns these two sequences
3. **samtools** converts `sam` file to `bam` file (*currently optional*)
4. **samtools** sorts `bam` file (*currently optional*)
5. **sam2pairwise** converts `sam` file to *pairwise* (`txt` file)


## Requirements

- **BWA**: http://bio-bwa.sourceforge.net
- **Samtools**: https://www.htslib.org
- **sam2pairwise**: https://github.com/mlafave/sam2pairwise

`sudo apt install bwa samtools`


## Run

###### Basic

- Navigate to `BWA/small/{...}` or `BWA/large{N}` folder
- Run `./run.sh`

###### Pairwise

`sam2pairwise < bwa_output.sam > bwa_output_pairwise.txt`


## Contents

- `large1/large_genome1.fasta`: [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/CP003305.1)
- `large1/large_genome2.fasta`: [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/CP000766.3)
---
- `large2/large_genome1.fasta`: [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP009625.1)
- `large2/large_genome2.fasta`: [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP007695.1)
---
- `large3/large_genome1.fasta`: [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP009626.1)
- `large3/large_genome2.fasta`: [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP007696.1)

## Usefull

https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform

https://ru.wikipedia.org/wiki/SAMtools

https://en.wikipedia.org/wiki/SAM_(file_format)

http://bio-bwa.sourceforge.net/bwa.shtml#3

http://bio-bwa.sourceforge.net/bwa.shtml#4

<p align="center">CIGAR</p>

![](src/CIGAR.png)