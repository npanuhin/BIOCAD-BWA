# import re
from string import ascii_uppercase
import matplotlib.pyplot as plt
# from random import randint
from Bio import SeqIO

FONT_SIZE = 9
DOT_SIZE = 1.5

BITWISE_FLAGS = [
    "template having multiple templates in sequencing (read is paired)",
    "each segment properly aligned according to the aligner (read mapped in proper pair)",
    "segment unmapped (read1 unmapped)",
    "next segment in the template unmapped (read2 unmapped)",
    "SEQ being reverse complemented (read1 reverse complemented)",
    "SEQ of the next segment in the template being reverse complemented (read2 reverse complemented)",
    "the first segment in the template (is read1)",
    "the last segment in the template (is read2)",
    "not primary alignment",
    "alignment fails quality checks",
    "PCR or optical duplicate",
    "supplementary alignment (e.g. aligner specific, could be a portion of a split read or a tied region)"
]


def prettifyNumber(num):
    return "{:,}".format(num)

def main(query_genome_path, ref_genome_path, sam_file_path):

    for record in SeqIO.parse(query_genome_path, "fasta"):
        query_genome_legnth = len(record.seq)
        break

    for record in SeqIO.parse(ref_genome_path, "fasta"):
        ref_genome_legnth = len(record.seq)
        break

    print("query_genome_legnth: ", query_genome_legnth)
    print("ref_genome_legnth: ", ref_genome_legnth)

    # sam_info = []
    sam_data = []
    query_genome_name, ref_genome_name = None, None

    with open(sam_file_path, 'r', encoding="utf-8") as sam_file:
        for line in sam_file:
            if line.startswith('@'):
                # sam_info.append(line)
                continue

            line = line.split()

            # Quality
            mapQ = int(line[4])
            quality = round((10 ** (mapQ / -10)) * 100, 6)

            # Start position
            position = int(line[3])

            # Flags
            flags_bit = int(line[1])
            flags = []
            for i in range(len(BITWISE_FLAGS) - 1, -1, -1):
                cur_flag, flags_bit = divmod(flags_bit, 2 ** i)
                if cur_flag:
                    flags.append(i)
            flags.sort()

            # Rif
            rid = line[9]
            rid_size = len(rid)
            if rid_size <= 10000:
                continue

            # CIGAR
            cigar = line[5]
            actions = []
            buff = ""
            for i in range(len(cigar)):
                if cigar[i] in ascii_uppercase:
                    actions.append([int(buff), cigar[i]])
                    # print("{}-{}|".format(buff, actions[i]), end="")
                    buff = ""
                else:
                    buff += cigar[i]

            # Names
            if ref_genome_name is None:
                ref_genome_name = line[0]
            if query_genome_name is None:
                query_genome_name = line[2]

            sam_data.append([position, flags, quality, rid_size, actions])


    sam_data.sort(key=lambda item: item[2])

    # fig, plot = plt.subplots()
    plt.rc("font", size=FONT_SIZE)               # controls default text sizes
    plt.rc("axes", titlesize=FONT_SIZE)          # fontsize of the axes title
    plt.rc("axes", labelsize=8)                  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=FONT_SIZE)         # fontsize of the tick labels
    plt.rc("ytick", labelsize=FONT_SIZE)         # fontsize of the tick labels
    plt.rc("legend", fontsize=FONT_SIZE)         # legend fontsize
    plt.rc("figure", titlesize=FONT_SIZE)        # fontsize of the figure title
    plt.ticklabel_format(style = "plain")
    plt.xticks(list(range(0, int(1e7), int(1e5))))
    plt.yticks(list(range(0, int(1e7), int(1e5))))
    plt.grid(which="major", linestyle='-', linewidth="0.5", color="red")
    # plt.minorticks_on()
    plt.grid(which="minor", linestyle=':', linewidth="0.5", color="black")
    plt.xlabel(query_genome_name)
    plt.ylabel(ref_genome_name)
    # plot.vlines(2, y.min(), y.max(), color = 'r')  # Вертикальные линии
    # plot.hlines(-4, -5, 5)  # Горизонтальные динии

    dots = []

    for position, flags, quality, rid_size, actions in sam_data:
        print("{} -> {} [Q:{}]".format(prettifyNumber(position), prettifyNumber(position + rid_size), quality))
        for flag in flags:
            if flag == 11:  # Currently flag №11 is disabled
                continue
            print(BITWISE_FLAGS[flag])

        # for length, action_type in actions:
        #     print("{}-{}|".format(length, action_type), end="")
        # print()

        cur_query_pos = 0
        cur_ref_pos = position
        cur_dots = []

        for length, action_type in actions:
            if action_type == 'S':
                cur_query_pos += length
                # cur_ref_pos += length

            elif action_type == 'H':
                cur_query_pos += length
                # cur_ref_pos += length

            elif action_type == 'M':
                for i in range(length):
                    cur_dots.append([cur_ref_pos, cur_query_pos])
                    cur_query_pos += 1
                    cur_ref_pos += 1

            elif action_type == 'I':
                for i in range(length):
                    cur_dots.append([cur_ref_pos, cur_query_pos])
                    # cur_query_pos += 1
                    cur_ref_pos += 1

            elif action_type == 'D':
                for i in range(length):
                    cur_dots.append([cur_ref_pos, cur_query_pos])
                    cur_query_pos += 1
                    # cur_ref_pos += 1

            else:
                raise "Unknown action type"


        if 4 in flags:
            for i in range(len(cur_dots)):
                cur_dots[i][1] = ref_genome_legnth - cur_dots[i][1]

        dots += cur_dots

        # plot.vlines(position, 0, 1e7, color = 'r')
        # plot.vlines(position + rid_size, 0, 1e7, color = 'g')

        # plot.fill_between([position, position + rid_size], [0, 0], [1e7, 1e7],
        #     alpha = 0.5,
        #     color = (randint(0, 255) / 255, randint(0, 255) / 255, randint(0, 255) / 255),
        #     linewidth = 1)

        print()

    dots = dots[::100]

    print("Dots count:", len(dots))

    plt.scatter(*zip(*dots), s=DOT_SIZE).set_label('')
    plt.show()


# ---SETTINGS--- #

query_genome_path = "samples/large3/large_genome1.fasta"
ref_genome_path = "samples/large3/large_genome2.fasta"
sam_file_path = "BWA/large3/bwa_output.sam"

# ---SETTINGS--- #

if __name__ == "__main__":
    main(query_genome_path, ref_genome_path, sam_file_path)