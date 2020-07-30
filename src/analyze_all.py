import os
from sys import path as sys_path
sys_path.append("../")

import sam_analyze

sam_analyze.GRID_SIZE = 1e2
sam_analyze.MIN_RID_SIZE = 1
sam_analyze.DOT_SKIP_RATE = 1
sam_analyze.DOT_SIZE = 0.1
sam_analyze.MIN_EVENT_SIZE = 0

os.chdir("../")

for foldername in os.listdir("tests/small"):

    sam_analyze.main(
        query_genome_path="samples/small/source.fasta",
        ref_genome_path="samples/small/{}.fasta".format(foldername),
        sam_file_path="BWA/small/{}/bwa_output.sam".format(foldername),
        output_to_file=True,
        output_folder="tests/small/{}".format(foldername)
    )


sam_analyze.DOT_SKIP_RATE = 5
sam_analyze.GRID_SIZE = 1e5
sam_analyze.MIN_RID_SIZE = 1e3
sam_analyze.MIN_EVENT_SIZE = 500
sam_analyze.DOT_SIZE = 0.01

for foldername in os.listdir("tests"):
    if foldername.strip("/").strip("\\") == "small":
        continue

    sam_analyze.main(
        query_genome_path="samples/{}/large_genome1.fasta".format(foldername),
        ref_genome_path="samples/{}/large_genome2.fasta".format(foldername),
        sam_file_path="BWA/{}/bwa_output.sam".format(foldername),
        output_to_file=True,
        output_folder="tests/{}".format(foldername)
    )
