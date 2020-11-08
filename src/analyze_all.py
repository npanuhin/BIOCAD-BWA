import os
from sys import path as sys_path

os.chdir("../")

sys_path.append(".")

import sam_analyze


def mkpath(*paths):
    return os.path.normpath(os.path.join(*paths))


# Small
SETTINGS = {
    "grid_size": 100,
    "min_rid_size": 1,
    "dot_skip_rate": 1,
    "dotsize": 0.1,
    "fontsize": 10,
    "figsize": (10, 7),

    "min_event_size": 3,
    "lines_join_size": 5,
    "line_min_size": 10
}

for foldername in os.listdir(mkpath("tests", "small")):
    break
    sam_analyze.analyze(
        query_genome_path="samples/small/source.fasta",
        ref_genome_path="samples/small/{}.fasta".format(foldername),
        sam_file_path="BWA/small/{}/bwa_output.sam".format(foldername),
        show_plot=False,
        output_folder="tests/small/{}".format(foldername),
        settings=SETTINGS.copy()
    )

# Large
SETTINGS = {
    "grid_size": int(1e5),
    "min_rid_size": int(1e3),
    "dot_skip_rate": 10,
    "dotsize": 0.1,
    "fontsize": 8,
    "figsize": (10, 7),

    "min_event_size": int(5e3),
    "lines_join_size": "$min_event_size + 3",
    "line_min_size": "$min_event_size"
}

for foldername in os.listdir("tests"):
    if not os.path.isdir(mkpath("tests", foldername)) or foldername.strip("/").strip("\\") == "small":
        continue

    # if foldername.strip("/").strip("\\") != "large06":
    #     continue

    sam_analyze.analyze(
        query_genome_path="samples/{}/large_genome1.fasta".format(foldername),
        ref_genome_path="samples/{}/large_genome2.fasta".format(foldername),
        sam_file_path="BWA/{}/bwa_output.sam".format(foldername),
        show_plot=False,
        output_folder="tests/{}".format(foldername),
        settings=SETTINGS.copy()
    )
