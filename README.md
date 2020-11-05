<h1 align="center">BWA and SAM analysis</h1>

> This repository is not intended to represent work, but rather to store and transmit data.

## Analysis status

<h3 align="center"><a href="https://github.com/npanuhin/BIOCAD_BWA/blob/master/tests/README.md">View current test results</a></h3>

✅ - Works as intended
⚠ - There are problems, but the solution is possible
❌ - There are problems that make the solution wrong

- ✅ large01
- ✅ large02
- ✅ large03
- ✅ large04
- ✅ large05
- ❌ large06
- ✅ large07
- 🔜 large08
- 🔜 large09
- ✅ small (BWA⚠)


## How BWA works now

1. **BWA** indexes two `fasta` sequences
2. **BWA** aligns these two sequences
3. **samtools** converts `sam` file to `bam` file (*currently disabled*)
4. **samtools** sorts `bam` file (*currently disabled*)
5. **sam2pairwise** converts `sam` file to *pairwise* (`txt` file) (*currently disabled*)

> For `SAM` and `pairwise` files *word wrap* should be disabled


#### Download software:

- **BWA**: http://bio-bwa.sourceforge.net
- **Samtools**: https://www.htslib.org
- **sam2pairwise**: https://github.com/mlafave/sam2pairwise

Or run `sudo apt install bwa samtools`


## Contents

1. [`large01/large_genome1.fasta`](./samples/large01 "Go to /samples/large01"): [Rickettsia rickettsii str. Brazil, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP003305.1)

   [`large01/large_genome2.fasta`](./samples/large01 "Go to /samples/large01"): [Rickettsia rickettsii str. Iowa, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP000766.3)
---
2. [`large02/large_genome1.fasta`](./samples/large02 "Go to /samples/large02"): [Brucella abortus 104M chromosome 1, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP009625.1)

   [`large02/large_genome2.fasta`](./samples/large02 "Go to /samples/large02"): [Brucella suis bv. 2 strain Bs143CITA chromosome I, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP007695.1)
---
3. [`large03/large_genome1.fasta`](./samples/large03 "Go to /samples/large03"): [Brucella abortus 104M chromosome 2, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP009626.1)

   [`large03/large_genome2.fasta`](./samples/large03 "Go to /samples/large03"): [Brucella suis bv. 2 strain Bs143CITA chromosome II, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP007696.1)
---
4. [`large04/large_genome1.fasta`](./samples/large04 "Go to /samples/large04"): [Brucella pinnipedialis B2/94 chromosome 2, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/CP002079)

   [`large04/large_genome2.fasta`](./samples/large04 "Go to /samples/large04"): [Brucella melitensis biovar Abortus 2308 chromosome II, complete sequence, strain 2308](https://www.ncbi.nlm.nih.gov/nuccore/AM040265.1)
---
5. [`large05/large_genome1.fasta`](./samples/large05 "Go to /samples/large05"): [Rickettsia rickettsii str. Iowa, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/864354655)

   [`large05/large_genome2.fasta`](./samples/large05 "Go to /samples/large05"): [Rickettsia prowazekii str. Madrid E, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/15603881)
---
6. [`large06/large_genome1.fasta`](./samples/large06 "Go to /samples/large06"): [Methanococcus maripaludis C5, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/134045046)

   [`large06/large_genome2.fasta`](./samples/large06 "Go to /samples/large06"): [Methanococcus maripaludis X1, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP002913.1)
---
7. [`large07/large_genome1.fasta`](./samples/large07 "Go to /samples/large07"): [Mycobacterium tuberculosis variant africanum GM041182, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_015758.1)

   [`large07/large_genome2.fasta`](./samples/large07 "Go to /samples/large07"): [Mycobacterium intracellulare ATCC 13950, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NC_016946.1)
---
8. [`large08/large_genome1.fasta`](./samples/large08 "Go to /samples/large08"): [Desulfurococcus kamchatkensis 1221n, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP001140.1)

   [`large08/large_genome2.fasta`](./samples/large08 "Go to /samples/large08"): [Desulfurococcus fermentans DSM 16532, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP003321.1)
---
9. [`large09/large_genome1.fasta`](./samples/large09 "Go to /samples/large09"): [Sulfolobus islandicus M.16.27, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP001401.1)

   [`large09/large_genome2.fasta`](./samples/large09 "Go to /samples/large09"): [Sulfolobus islandicus REY15A, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP002425.1)