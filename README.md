<h1 align="center">BWA and SAM analysis</h1>

> This repository is not intended to represent work, but rather to store and transmit data.

## Analysis status

<h3 align="center"><a href="https://github.com/npanuhin/BIOCAD_BWA/blob/master/tests/README.md">View current test results</a></h3>

✅ - Works as intended
⚠ - There are problems, but the solution is possible
❌ - There are problems that make the solution wrong

- ✅ large1
- ✅ large2
- ✅ large3
- ✅ large4
- ✅ large5
- ❌ large6
- ✅ large7
- 🔜 large8
- 🔜 large9
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

1. [`large1/large_genome1.fasta`](./samples/large1 "Go to /samples/large1"): [Rickettsia rickettsii str. Brazil, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP003305.1)

   [`large1/large_genome2.fasta`](./samples/large1 "Go to /samples/large1"): [Rickettsia rickettsii str. Iowa, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP000766.3)
---
2. [`large2/large_genome1.fasta`](./samples/large2 "Go to /samples/large2"): [Brucella abortus 104M chromosome 1, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP009625.1)

   [`large2/large_genome2.fasta`](./samples/large2 "Go to /samples/large2"): [Brucella suis bv. 2 strain Bs143CITA chromosome I, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP007695.1)
---
3. [`large3/large_genome1.fasta`](./samples/large3 "Go to /samples/large3"): [Brucella abortus 104M chromosome 2, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP009626.1)

   [`large3/large_genome2.fasta`](./samples/large3 "Go to /samples/large3"): [Brucella suis bv. 2 strain Bs143CITA chromosome II, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP007696.1)
---
4. [`large4/large_genome1.fasta`](./samples/large4 "Go to /samples/large4"): [Brucella pinnipedialis B2/94 chromosome 2, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/CP002079)

   [`large4/large_genome2.fasta`](./samples/large4 "Go to /samples/large4"): [Brucella melitensis biovar Abortus 2308 chromosome II, complete sequence, strain 2308](https://www.ncbi.nlm.nih.gov/nuccore/AM040265.1)
---
5. [`large5/large_genome1.fasta`](./samples/large5 "Go to /samples/large5"): [Rickettsia rickettsii str. Iowa, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/864354655)

   [`large5/large_genome2.fasta`](./samples/large5 "Go to /samples/large5"): [Rickettsia prowazekii str. Madrid E, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/15603881)
---
6. [`large6/large_genome1.fasta`](./samples/large6 "Go to /samples/large6"): [Methanococcus maripaludis C5, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/134045046)

   [`large6/large_genome2.fasta`](./samples/large6 "Go to /samples/large6"): [Methanococcus maripaludis X1, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP002913.1)
---
7. [`large7/large_genome1.fasta`](./samples/large7 "Go to /samples/large7"): [Mycobacterium tuberculosis variant africanum GM041182, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_015758.1)

   [`large7/large_genome2.fasta`](./samples/large7 "Go to /samples/large7"): [Mycobacterium intracellulare ATCC 13950, complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/NC_016946.1)
---
8. [`large8/large_genome1.fasta`](./samples/large8 "Go to /samples/large8"): [Desulfurococcus kamchatkensis 1221n, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP001140.1)

   [`large8/large_genome2.fasta`](./samples/large8 "Go to /samples/large8"): [Desulfurococcus fermentans DSM 16532, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP003321.1)
---
9. [`large9/large_genome1.fasta`](./samples/large9 "Go to /samples/large9"): [Sulfolobus islandicus M.16.27, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP001401.1)

   [`large9/large_genome2.fasta`](./samples/large9 "Go to /samples/large9"): [Sulfolobus islandicus REY15A, complete genome](https://www.ncbi.nlm.nih.gov/nuccore/CP002425.1)