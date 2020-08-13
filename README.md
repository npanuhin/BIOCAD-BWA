<h1 align="center">BWA</h1>

> This repository is not intended to represent work, but rather to store and transmit data.

For `SAM` and `pairwise` files *word wrap* should be disabled


## How BWA works now

1. **BWA** indexes two `fasta` sequences
2. **BWA** aligns these two sequences
3. **samtools** converts `sam` file to `bam` file (*currently optional*)
4. **samtools** sorts `bam` file (*currently optional*)
5. **sam2pairwise** converts `sam` file to *pairwise* (`txt` file) (*currently optional*)


#### Where to download software?

- **BWA**: http://bio-bwa.sourceforge.net
- **Samtools**: https://www.htslib.org
- **sam2pairwise**: https://github.com/mlafave/sam2pairwise

`sudo apt install bwa samtools`


## Contents

- `large1/large_genome1.fasta`: [Rickettsia rickettsii str. Brazil, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP003305.1)
- `large1/large_genome2.fasta`: [Rickettsia rickettsii str. Iowa, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP000766.3)
---
- `large2/large_genome1.fasta`: [Brucella abortus 104M chromosome 1, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP009625.1)
- `large2/large_genome2.fasta`: [Brucella suis bv. 2 strain Bs143CITA chromosome I, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP007695.1)
---
- `large3/large_genome1.fasta`: [Brucella abortus 104M chromosome 2, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP009626.1)
- `large3/large_genome2.fasta`: [Brucella suis bv. 2 strain Bs143CITA chromosome II, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP007696.1)
---
- `large4/large_genome1.fasta`: [Brucella pinnipedialis B2/94 chromosome 2, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/CP002079)
- `large4/large_genome2.fasta`: [Brucella melitensis biovar Abortus 2308 chromosome II, complete sequence, strain 2308](https://www.ncbi.nlm.nih.gov/nuccore/AM040265.1)
---
- `large5/large_genome1.fasta`: [Rickettsia rickettsii str. Iowa, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/864354655)
- `large5/large_genome2.fasta`: [Rickettsia prowazekii str. Madrid E, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/15603881)
---
- `large6/large_genome1.fasta`: [Methanococcus maripaludis C5, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/134045046)
- `large6/large_genome2.fasta`: [Methanococcus maripaludis X1, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP002913.1)
---
- `large7/large_genome1.fasta`: [Mycobacterium tuberculosis variant africanum GM041182, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_015758.1)
- `large7/large_genome2.fasta`: [Mycobacterium intracellulare ATCC 13950, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NC_016946.1)