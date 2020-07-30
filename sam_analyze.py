from string import ascii_uppercase
import matplotlib.pyplot as plt
from Bio import SeqIO
# from matplotlib.patches import ConnectionStyle
import matplotlib.lines as mlines
import os
# from random import randint
# import re

INT_MAX = int(1e9) + 7


# !!! X - query, Y - ref !!!

# Small
# GRID_SIZE = 1e2
# MIN_RID_SIZE = 1e2
# DOT_SKIP_RATE = 1
# DOT_SIZE = 0.1
# MIN_EVENT_SIZE = 0

# Large
GRID_SIZE = 1e5
MIN_RID_SIZE = 1e3
DOT_SKIP_RATE = 10
DOT_SIZE = 0.1
MIN_EVENT_SIZE = 500

EPSILON = 1e4

FONT_SIZE = 10
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

# ---SETTINGS--- #

query_genome_path = "samples/large4/large_genome1.fasta"
ref_genome_path = "samples/large4/large_genome2.fasta"
sam_file_path = "BWA/large4/bwa_output.sam"

output_to_file = True

output_folder = "tests/large4"

# ---SETTINGS--- #


def mkpath(*paths):
    return os.path.normpath(os.path.join(*paths))


def prettifyNumber(num):
    return "{:,}".format(num)


def drawline(fig, p1, p2):
    fig.add_artist(mlines.Line2D(p1, p2))


def equalE(value1, value2):
    return value1 - EPSILON < value2 < value1 + EPSILON


# Plot setup
def setPlot(nameX, nameY):
    plt.rc("font", size=FONT_SIZE)               # controls default text sizes
    plt.rc("axes", titlesize=FONT_SIZE)          # fontsize of the axes title
    plt.rc("axes", labelsize=8)                  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=FONT_SIZE)         # fontsize of the tick labels
    plt.rc("ytick", labelsize=FONT_SIZE)         # fontsize of the tick labels
    plt.rc("legend", fontsize=FONT_SIZE)         # legend fontsize
    plt.rc("figure", titlesize=FONT_SIZE)        # fontsize of the figure title
    plt.ticklabel_format(style="plain")
    plt.xticks(list(range(0, int(GRID_SIZE * 100), int(GRID_SIZE))))
    plt.yticks(list(range(0, int(GRID_SIZE * 100), int(GRID_SIZE))))
    plt.xticks(rotation=30)
    plt.grid(which="major", linestyle='-', linewidth="1", alpha=0.1, color="black")
    # plt.minorticks_on()
    # plt.grid(which="minor", linestyle=':', linewidth="1", color="black")
    plt.xlabel(nameX)
    plt.ylabel(nameY)
    # plot.vlines(2, y.min(), y.max(), color = 'r')  # Вертикальные линии
    # plot.hlines(-4, -5, 5)  # Горизонтальные динии


def main(query_genome_path, ref_genome_path, sam_file_path, output_to_file, output_folder):

    for record in SeqIO.parse(query_genome_path, "fasta"):
        query_genome_length = len(record.seq)
        break

    for record in SeqIO.parse(ref_genome_path, "fasta"):
        ref_genome_length = len(record.seq)
        break

    print("query_genome_length: ", query_genome_length)
    print("ref_genome_length: ", ref_genome_length)
    print()

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
            if rid_size <= MIN_RID_SIZE:
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

    plt.figure("Main", figsize=(16, 11.2))
    setPlot(query_genome_name, ref_genome_name)

    # connectionstyle_I = ConnectionStyle("Angle3", angleA=90, angleB=0)
    # connectionstyle_D = ConnectionStyle("Angle3", angleA=0, angleB=90)

    all_actions = []
    large_actions = []

    for position, flags, quality, rid_size, actions in sam_data:
        print("{} -> {} [Q:{}]".format(prettifyNumber(position), prettifyNumber(position + rid_size), quality))
        for flag in flags:
            if flag == 11:  # Currently flag №11 is disabled
                continue
            print(BITWISE_FLAGS[flag])

        # for length, action_type in actions:
        #     print("{}-{}|".format(length, action_type), end="")
        # print()

        cur_query_pos = position
        cur_ref_pos = 0
        start_query_pos, start_ref_pos, end_query_pos, end_ref_pos = INT_MAX, INT_MAX, INT_MAX, INT_MAX

        for length, action_type in actions:

            if action_type not in ('S', 'H'):
                start_query_pos = min(start_query_pos, cur_query_pos)
                start_ref_pos = min(start_ref_pos, cur_ref_pos)

            if action_type == 'S':
                cur_ref_pos += length

            elif action_type == 'H':
                cur_ref_pos += length

            elif action_type == 'M':
                all_actions.append([cur_query_pos, cur_ref_pos, length, action_type, flags])
                cur_query_pos += length
                cur_ref_pos += length

            elif action_type == 'I':
                all_actions.append([cur_query_pos, cur_ref_pos, length, action_type, flags])
                cur_ref_pos += length

            elif action_type == 'D':
                all_actions.append([cur_query_pos, cur_ref_pos, length, action_type, flags])
                cur_query_pos += length

            else:
                raise "Unknown action type"

            # if action_type not in ('S', 'H'):
            #     end_ref_pos = min(end_ref_pos, cur_query_pos)
            #     end_query_pos = min(end_query_pos, cur_ref_pos)

        if 4 in flags:
            length = sum([action[0] for action in actions if action[1] not in ('S', 'H')])

            if length >= MIN_EVENT_SIZE:
                large_actions.append(['R', start_query_pos, start_ref_pos, length, flags])

        print()

    all_actions.sort(key=lambda action: action[0])
    last_query_end, last_ref_end, last_length, last_flags = None, None, None, None
    dots, ghost_dots = [], []

    for action_index in range(len(all_actions)):

        cur_query_pos, cur_ref_pos, length, action_type, flags = all_actions[action_index]

        if action_index > 0:

            ref_length = abs(cur_ref_pos - last_ref_end)
            query_length = abs(cur_query_pos - last_query_end)

            if ref_length > query_length:

                if ref_length >= MIN_EVENT_SIZE:
                    large_actions.append(['D', last_query_end, last_ref_end, ref_length])

                plt.plot([last_query_end, cur_query_pos], [last_ref_end, cur_ref_pos], color="red")

            else:

                if query_length >= MIN_EVENT_SIZE:
                    large_actions.append(['I', last_query_end, last_ref_end, ref_length])

                plt.plot([last_query_end, cur_query_pos], [last_ref_end, cur_ref_pos], color="green")

                # plt.axvline(x=cur_ref_pos, color="#ff0")
                # if length >= MIN_EVENT_SIZE:
                #     plt.plot([cur_ref_pos, cur_ref_pos + 50], [cur_query_pos, cur_query_pos], color="green")
                #     plt.plot([cur_ref_pos, cur_ref_pos + 50], [cur_query_pos + length, cur_query_pos + length], color="green")
                #     plt.plot([cur_ref_pos + 50, cur_ref_pos + 50], [cur_query_pos + length, cur_query_pos], color="green")

                #     plt.annotate("insertion", xy=(cur_ref_pos + 50, cur_query_pos + length / 2),
                #                  arrowprops=dict(arrowstyle="-", connectionstyle=connectionstyle_I, color="green"),
                #                  xycoords='data',
                #                  xytext=(10, -20),
                #                  color="green",
                #                  textcoords='offset points'
                #                  )

        if action_type == 'M':
            for i in range(length):
                if 4 in flags:
                    dots.append([cur_query_pos, ref_genome_length - cur_ref_pos])
                    ghost_dots.append([cur_query_pos, cur_ref_pos])
                else:
                    dots.append([cur_query_pos, cur_ref_pos])
                cur_query_pos += 1
                cur_ref_pos += 1

        elif action_type == 'I':
            plt.plot([cur_query_pos, cur_query_pos], [cur_ref_pos, cur_ref_pos + length], color="green")

            if length >= MIN_EVENT_SIZE:
                large_actions.append([action_type, cur_query_pos, cur_ref_pos, length])
            #     plt.plot([cur_ref_pos, cur_ref_pos + 50], [cur_query_pos + length, cur_query_pos + length], color="green")
            #     plt.plot([cur_ref_pos + 50, cur_ref_pos + 50], [cur_query_pos + length, cur_query_pos], color="green")
            #     plt.annotate("insertion", xy=(cur_ref_pos + 50, cur_query_pos + length / 2),
            #                  arrowprops=dict(arrowstyle="-", connectionstyle=connectionstyle_I, color="green"),
            #                  xycoords='data',
            #                  xytext=(10, -20),
            #                  color="green",
            #                  textcoords='offset points'
            #                  )

            cur_ref_pos += length

        elif action_type == 'D':
            plt.plot([cur_query_pos, cur_query_pos + length], [cur_ref_pos, cur_ref_pos], color="red")

            if length >= MIN_EVENT_SIZE:
                large_actions.append([action_type, cur_query_pos, cur_ref_pos, length])
            #     plt.plot([cur_ref_pos, cur_ref_pos], [cur_query_pos, cur_query_pos - 50], color="red")
            #     plt.plot([cur_ref_pos + length, cur_ref_pos + length], [cur_query_pos, cur_query_pos - 50], color="red")
            #     plt.plot([cur_ref_pos + length, cur_ref_pos], [cur_query_pos - 50, cur_query_pos - 50], color="red")
            #     plt.annotate("deletion", xy=(cur_ref_pos + length / 2, cur_query_pos - 50),
            #                  arrowprops=dict(arrowstyle="-", connectionstyle=connectionstyle_D, color="red"),
            #                  xycoords='data',
            #                  xytext=(10, -20),
            #                  color="red",
            #                  textcoords='offset points'
            #                  )

            cur_query_pos += length

        last_ref_end = cur_ref_pos
        last_query_end = cur_query_pos
        # last_length = length
        # last_flags = flags

    dots = dots[::DOT_SKIP_RATE]
    ghost_dots = ghost_dots[::DOT_SKIP_RATE]

    print("Dots count: {} + {}".format(prettifyNumber(len(dots)), prettifyNumber(len(ghost_dots))))

    plt.scatter(*zip(*dots), s=DOT_SIZE)
    if ghost_dots:
        plt.scatter(*zip(*ghost_dots), s=DOT_SIZE, color="#ccc")
    plt.tight_layout()

    if output_to_file:
        print("Saving plot...")
        plt.savefig(mkpath(output_folder, "sam_analyze.png"), dpi=300)
    else:
        print("Showing plot...")
        plt.show()

    plt.cla()
    plt.clf()
    plt.close()

# ------------------------------------------------------------------------------------------------

    large_actions.sort(key=lambda action: action[3])

    if not os.path.exists(mkpath(output_folder, "history")):
        os.mkdir(mkpath(output_folder, "history"))

    for filename in os.listdir(mkpath(output_folder, "history")):
        os.remove(mkpath(output_folder, "history", filename))

    return  # TODO

    action_index = 1
    for action_index in range(len(large_actions)):
        action = large_actions[action_index]
        action_type, start_ref_pos, start_query_pos, length = action[0], action[1], action[2], action[3]

        plt.figure("large_action{}".format(large_actions), figsize=(16, 11.2))
        setPlot(query_genome_name, ref_genome_name)

        new_dots = []

        if action_type == 'R':
            for dot_x, dot_y in dots:
                if start_ref_pos <= dot_x <= start_ref_pos + length:
                    dot_y

        dots = new_dots

        # TODO

        # plt.scatter(*zip(*dots), s=DOT_SIZE)
        plt.tight_layout()
        print("Saving large action # {}...".format(action_index))
        plt.savefig(mkpath(output_folder, "history", str(action_index) + ".png"), dpi=300)

        plt.cla()
        plt.clf()
        plt.close()


if __name__ == "__main__":
    main(query_genome_path, ref_genome_path, sam_file_path, output_to_file, output_folder)
