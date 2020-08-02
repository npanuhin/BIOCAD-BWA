from string import ascii_uppercase
import matplotlib.pyplot as plt
from Bio.SeqIO.FastaIO import SimpleFastaParser
# from matplotlib.patches import ConnectionStyle
# import matplotlib.lines as mlines
import os
# from collections import deque
# from threading import Thread
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
MIN_EVENT_SIZE = 1e3

ROTATION_JOIN_SIZE = 1e4

FIGSIZE = (10, 7)
FONT_SIZE = 10

CIGAR_FLAGS = [
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

query_genome_path = "samples/large3/large_genome1.fasta"
ref_genome_path = "samples/large3/large_genome2.fasta"
sam_file_path = "BWA/large3/bwa_output.sam"

show_plot = True

output_folder = "tests/large3"

# ---SETTINGS--- #


def mkpath(*paths):
    return os.path.normpath(os.path.join(*paths))


def prettifyNumber(num):
    return "{:,}".format(num)


def equalE(value1, value2, E):
    '''Epsilon comparison'''
    return value1 - E < value2 < value1 + E


class Plot:
    def __init__(self, title=None, nameX=None, nameY=None, figsize=FIGSIZE):
        self.fig = plt.figure(title, figsize)
        self.ax = self.fig.add_subplot()

        # self.fig.rc("font", size=FONT_SIZE)               # controls default text sizes
        # self.fig.rc("axes", titlesize=FONT_SIZE)          # fontsize of the axes title
        # self.ax.rc("axes", labelsize=8)                  # fontsize of the x and y labels
        # self.ax.rc("xtick", labelsize=FONT_SIZE)         # fontsize of the tick labels
        # self.ax.rc("ytick", labelsize=FONT_SIZE)         # fontsize of the tick labels
        # self.ax.rc("legend", fontsize=FONT_SIZE)         # legend fontsize
        # self.ax.rc("figure", titlesize=FONT_SIZE)        # fontsize of the figure title
        self.ax.ticklabel_format(style="plain")
        self.ax.set_xticks(list(range(0, int(GRID_SIZE * 1000), int(GRID_SIZE))))
        self.ax.set_yticks(list(range(0, int(GRID_SIZE * 1000), int(GRID_SIZE))))
        self.ax.tick_params(axis="x", which="major", labelsize=FONT_SIZE, rotation=30)
        self.ax.tick_params(axis="y", which="major", labelsize=FONT_SIZE)
        self.ax.grid(which="major", linestyle='-', linewidth="1", alpha=0.1, color="black")

        # self.ax.minorticks_on()
        # self.ax.grid(which="minor", linestyle=':', linewidth="1", color="black")

        if nameX is not None:
            self.ax.set_xlabel(nameX)
        if nameY is not None:
            self.ax.set_ylabel(nameY)

    def __del__(self):
        self.close()

    def scatter(self, dots, s=DOT_SIZE, *args, **kwargs):
        self.ax.scatter(*zip(*dots), s=DOT_SIZE, *args, **kwargs)

    def line(self, x1, y1, x2, y2, *args, **kwargs):
        self.ax.plot([x1, x2], [y1, y2], *args, **kwargs)

    def save(self, path):
        self.fig.savefig(path, dpi=100)

    def close(self):
        plt.close(self.fig)


def main(query_genome_path, ref_genome_path, sam_file_path, show_plot, output_folder):

    with open(query_genome_path, 'r', encoding="utf-8") as file:
        for _, sequence in SimpleFastaParser(file):
            query_genome_length = len(sequence)
            break

    with open(ref_genome_path, 'r', encoding="utf-8") as file:
        for _, sequence in SimpleFastaParser(file):
            ref_genome_length = len(sequence)
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
            # mapQ = int(line[4])
            # quality = round((10 ** (mapQ / -10)) * 100, 6)

            # Start position
            position = int(line[3])

            # Flags
            flags_bit = int(line[1])
            flags = []
            for i in range(len(CIGAR_FLAGS) - 1, -1, -1):
                cur_flag, flags_bit = divmod(flags_bit, 2 ** i)
                if cur_flag:
                    flags.append(i)
            flags.sort()

            # Rid
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

            sam_data.append([position, flags, rid_size, actions])

    # sam_data.sort(key=lambda item: item[2])

    # return
# ========================================================================================================================================

    all_actions = []

    for position, flags, rid_size, actions in sam_data:
        print("{} -> {}".format(prettifyNumber(position), prettifyNumber(position + rid_size)))
        for flag in flags:
            # if flag == 11:  # Currently flag №11 is disabled for showing
            #     continue
            print(CIGAR_FLAGS[flag])

        # for length, action_type in actions:
        #     print("{}-{}|".format(length, action_type), end="")
        # print()

        cur_query_pos = position
        cur_ref_pos = 0
        start_query_pos, start_ref_pos, end_query_pos, end_ref_pos = INT_MAX, INT_MAX, 0, 0

        for length, action_type in actions:

            if action_type not in ('S', 'H'):
                start_query_pos = min(start_query_pos, cur_query_pos)
                start_ref_pos = min(start_ref_pos, cur_ref_pos)
                end_query_pos = max(end_query_pos, cur_query_pos)
                end_ref_pos = max(end_ref_pos, cur_ref_pos)

            if action_type == 'S':
                cur_ref_pos += length

            elif action_type == 'H':
                cur_ref_pos += length

            elif action_type == 'M':
                cur_query_pos += length
                cur_ref_pos += length

            elif action_type == 'I':
                cur_ref_pos += length

            elif action_type == 'D':
                cur_query_pos += length

            else:
                raise "Unknown action type"

        rotation_center = None
        if 4 in flags:
            rotation_center = (start_ref_pos + end_ref_pos) / 2

        cur_query_pos = position
        cur_ref_pos = 0

        for length, action_type in actions:

            if action_type == 'S':
                cur_ref_pos += length

            elif action_type == 'H':
                cur_ref_pos += length

            elif action_type == 'M':
                all_actions.append([cur_query_pos, cur_ref_pos, length, action_type, rotation_center])
                cur_query_pos += length
                cur_ref_pos += length

            elif action_type == 'I':
                all_actions.append([cur_query_pos, cur_ref_pos, length, action_type, rotation_center])
                cur_ref_pos += length

            elif action_type == 'D':
                all_actions.append([cur_query_pos, cur_ref_pos, length, action_type, rotation_center])
                cur_query_pos += length

        print()

    # return
# ========================================================================================================================================
    # Join rotations for all_actions + create rotations for large_actions

    rotations = []

    def getActionRefEnd(action):
        start_query_pos, start_ref_pos, length, action_type, rotation_center = action
        if action_type == 'M':
            return start_ref_pos + length
        elif action_type == 'I':
            return start_ref_pos
        elif action_type == 'D':
            return start_ref_pos + length

    def getActionQueryEnd(action):
        start_query_pos, start_ref_pos, length, action_type, rotation_center = action
        if action_type == 'M':
            return start_query_pos + length
        elif action_type == 'I':
            return start_query_pos + length
        elif action_type == 'D':
            return start_query_pos

    all_actions.sort(key=lambda action: action[0])

    rotation_block_start, rotation_block_end = 0, 0
    last_rotation_ref_end = None
    for i in range(len(all_actions)):
        action = all_actions[i]
        start_query_pos, start_ref_pos, length, action_type, rotation_center = action
        if rotation_center is None:
            continue

        if last_rotation_ref_end is None:
            rotation_block_start = i

        elif last_rotation_ref_end + ROTATION_JOIN_SIZE >= start_query_pos:
            pass

        else:
            new_rotation_center = (all_actions[rotation_block_start][1] + last_rotation_ref_end) / 2
            for j in range(rotation_block_start, rotation_block_end + 1):
                all_actions[j][4] = new_rotation_center

            block_length = getActionQueryEnd(all_actions[rotation_block_end]) - all_actions[rotation_block_start][0]
            if block_length >= MIN_EVENT_SIZE:
                rotations.append(['R', all_actions[rotation_block_start][0], all_actions[rotation_block_start][1], block_length, ref_genome_length - new_rotation_center])

            rotation_block_start = i

        rotation_block_end = i
        last_rotation_ref_end = getActionRefEnd(action)

    if last_rotation_ref_end is not None:
        new_rotation_center = (all_actions[rotation_block_start][1] + last_rotation_ref_end) / 2
        for j in range(rotation_block_start, rotation_block_end + 1):
            all_actions[j][4] = new_rotation_center

    block_length = getActionQueryEnd(all_actions[rotation_block_end]) - all_actions[rotation_block_start][0]
    if block_length >= MIN_EVENT_SIZE:
        rotations.append(['R', all_actions[rotation_block_start][0], all_actions[rotation_block_start][1], block_length, ref_genome_length - new_rotation_center])

    # return
    # ========================================================================================================================================

    all_actions.sort(key=lambda action: action[0])

    # if len(all_actions) > 1:

    #     all_actions_circle = deque(all_actions)

    #     while True:
    #         first_action = all_actions_circle[0]
    #         first_action_type, first_action_ref_pos, first_action_length = first_action[3], first_action[1], first_action[2]
    #         last_action = all_actions_circle[-1]
    #         last_action_type, last_action_ref_pos, last_action_length = last_action[3], last_action[1], last_action[2]

    #         first_ref_pos = first_action_ref_pos
    #         last_ref_pos = last_action_ref_pos

    #         if last_action_type == 'M':
    #             last_ref_pos += last_action_length
    #         elif last_action_type == 'I':
    #             last_ref_pos += last_action_length
    #         elif last_action_type == 'D':
    #             pass

    #         first_action_X_length = 0 if first_action_type == 'I' else first_action_length
    #         last_action_X_length = 0 if last_action_type == 'I' else last_action_length

    #         print(first_ref_pos, last_ref_pos)

    #         if equalE(first_ref_pos, last_ref_pos, 1e4):
    #             all_actions_circle.append(all_actions_circle.popleft())
    #             for i in range(len(all_actions) - 1):
    #                 all_actions_circle[i][0] -= first_action_X_length
    #             all_actions_circle[-1][0] = all_actions_circle[-2][0] + last_action_X_length
    #         else:
    #             break

    #     all_actions = list(all_actions_circle)

    # return
# ========================================================================================================================================

    plot = Plot("Main", query_genome_name, ref_genome_name)
    last_query_end, last_ref_end, last_rotated_ref_end = None, None, None
    dots, ghost_dots, rotated_dots, large_actions = [], [], [], []

    for action_index in range(len(all_actions)):

        start_query_pos, start_ref_pos, length, action_type, rotation_center = all_actions[action_index]
        cur_query_pos, cur_ref_pos = start_query_pos, start_ref_pos

        if rotation_center is None:
            rotated = lambda cur_ref_pos: cur_ref_pos
        else:
            rotated = lambda cur_ref_pos: (ref_genome_length - cur_ref_pos) + ((ref_genome_length - rotation_center) - (ref_genome_length - cur_ref_pos)) * 2

        if action_index > 0:
            ref_gap = last_rotated_ref_end - rotated(cur_ref_pos)
            # print(last_rotated_ref_end, rotated(cur_ref_pos), ref_gap)
            query_gap = cur_query_pos - last_query_end

            # The opposite gap !!!

            if abs(query_gap) + abs(ref_gap) > 0:

                color_pos = abs(ref_gap) / (abs(query_gap) + abs(ref_gap))

                plot.line(last_query_end, last_rotated_ref_end, cur_query_pos, rotated(cur_ref_pos), color=(1 - color_pos, color_pos, 0))

            # if query_gap > ref_gap and query_gap > 1:

            #     # if ref_gap >= MIN_EVENT_SIZE:
            #     #     large_actions.append(['D', last_query_end, last_ref_end, cur_ref_pos - last_ref_end])

            #     plot.line(last_query_end, last_ref_end, cur_query_pos, cur_ref_pos, color="#f00")

            # elif ref_gap > 1:

            #     # if query_gap >= MIN_EVENT_SIZE:
            #     #     large_actions.append(['I', last_query_end, last_ref_end, cur_query_pos - last_query_end])

            #     plot.line(last_query_end, last_ref_end, cur_query_pos, cur_ref_pos, color="#0f0")

            if abs(query_gap) >= MIN_EVENT_SIZE:  # The opposite gap !!!
                large_actions.append(['D', last_query_end, last_ref_end, query_gap])

            if abs(ref_gap) >= MIN_EVENT_SIZE:  # The opposite gap !!!
                # if last_query_end == 1098909:
                #     print("Here")
                print(last_query_end, last_ref_end, last_rotated_ref_end, rotated(cur_ref_pos), ref_gap)
                large_actions.append(['I', last_query_end, last_ref_end, ref_gap])

        if action_type == 'M':
            for i in range(length):
                if rotation_center is None:
                    dots.append([cur_query_pos, cur_ref_pos])
                    rotated_dots.append([cur_query_pos, cur_ref_pos])
                else:
                    dots.append([cur_query_pos, ref_genome_length - cur_ref_pos])
                    rotated_dots.append([cur_query_pos, rotated(cur_ref_pos)])
                    ghost_dots.append([cur_query_pos, rotated(cur_ref_pos)])

                cur_query_pos += 1
                cur_ref_pos += 1

        elif action_type == 'I':

            plot.line(cur_query_pos, rotated(cur_ref_pos), cur_query_pos, rotated(cur_ref_pos) + length, color="#0f0")

            if length >= MIN_EVENT_SIZE:
                large_actions.append([action_type, cur_query_pos, rotated(cur_ref_pos), length])

            cur_ref_pos += length

        elif action_type == 'D':
            plot.line(cur_query_pos, rotated(cur_ref_pos), cur_query_pos + length, rotated(cur_ref_pos), color="#f00")

            if length >= MIN_EVENT_SIZE:
                large_actions.append([action_type, cur_query_pos, rotated(cur_ref_pos), length])

            cur_query_pos += length

        if last_query_end is None or cur_query_pos >= last_query_end:
            last_ref_end = cur_ref_pos
            last_query_end = cur_query_pos
            last_rotated_ref_end = rotated(cur_ref_pos)

    dots = dots[::DOT_SKIP_RATE]
    ghost_dots = ghost_dots[::DOT_SKIP_RATE]
    rotated_dots = rotated_dots[::DOT_SKIP_RATE]

    print("Dots count: {} + {}".format(prettifyNumber(len(dots)), prettifyNumber(len(ghost_dots))))

    plot.scatter(dots)
    if ghost_dots:
        plot.scatter(ghost_dots, color="#ccc")
    plot.fig.tight_layout()

    print("Saving plot...")
    plot.save(mkpath(output_folder, "sam_analyze.png"))

    if show_plot:
        print("Showing plot...")
        plt.show()

    del plot

    # return
# ========================================================================================================================================
    # Join rotations for large_actions

    # rotations = sorted([action for action in large_actions if action[0] == 'R'], key=lambda action: action[1])

    # for i in range(len(rotations) - 1, 0, -1):
    #     last_query_pos, last_ref_pos, last_length = rotations[i - 1][1:4]
    #     query_pos, ref_pos, length = rotations[i][1:4]

    #     print(rotations[i])

    #     if last_query_pos + last_length + ROTATION_JOIN_SIZE >= query_pos:
    #         rotations[i - 1][3] = (query_pos + length) - last_query_pos  # last_length
    #         del rotations[i]

    # rotations.sort(key=lambda action: action[1])

    # print()
    # for rotation in rotations:
    #     print(rotation)
    # print()

    # large_actions = rotations + [action for action in large_actions if action[0] != 'R']

    dot_rotation_centers = [None] * len(dots)

    for action_type, start_query_pos, start_ref_pos, length, rotation_center in rotations:
        for i in range(len(dots)):
            if start_query_pos <= dots[i][0] <= start_query_pos + length:
                dot_rotation_centers[i] = rotation_center

    large_actions = rotations + large_actions

    large_actions.sort(key=lambda action: -abs(action[3]))

    # large_actions = [action for action in large_actions if action[0] == 'R'] + \
    #                 [action for action in large_actions if action[0] != 'R']

    print(large_actions)

    # return
# ========================================================================================================================================

    if not os.path.exists(mkpath(output_folder, "history")):
        os.mkdir(mkpath(output_folder, "history"))

    for filename in os.listdir(mkpath(output_folder, "history")):
        os.remove(mkpath(output_folder, "history", filename))

    large_actions = [['P', 0, 0, 0]] + large_actions

    for action_index in range(len(large_actions)):
        action = large_actions[action_index]
        print(action)
        # continue
        action_type, start_query_pos, start_ref_pos, length = action[0:4]

        # length - направленная длина!!!

        action_plot = Plot("large_action{}".format(action_index + 1), nameX=query_genome_name, nameY=ref_genome_name)

        if action_type == 'R':

            rotated = lambda cur_ref_pos, i: (ref_genome_length - cur_ref_pos) + (dot_rotation_centers[i] - (ref_genome_length - cur_ref_pos)) * 2
            # rotated = lambda cur_ref_pos: ref_genome_length - (action[4] + (action[4] - cur_ref_pos))

            for i in range(len(dots)):
                if start_query_pos <= dots[i][0] <= start_query_pos + length:
                    dots[i] = rotated_dots[i].copy()
                    # dots[i][1] = ref_genome_length - (dot_rotation_centers[i] + (dot_rotation_centers[i] - dots[i][1]))

        elif action_type == 'D':
            for i in range(len(dots)):
                if dots[i][0] >= start_query_pos + abs(length):
                    dots[i][0] -= length
            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos + abs(length):
                    large_actions[i][1] -= abs(length)

        elif action_type == 'I':
            print("I: query[{}] ref[{}] length[{}]".format(start_query_pos, start_ref_pos, length))
            for i in range(len(dots)):
                if dots[i][0] >= start_query_pos:
                    dots[i][1] += length
                    rotated_dots[i][1] += length

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] += length

                # if large_actions[i][0] == 'R':
                #     print("Was:", large_actions[i][4], "New:", large_actions[i][4] - length)
                #     large_actions[i][4] -= length

            bottom = INT_MAX
            for dot_x, dot_y in dots:
                bottom = min(bottom, dot_y)

            for i in range(len(dots)):
                dots[i][1] -= bottom
                rotated_dots[i][1] -= bottom

            # print(bottom, length)

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] -= bottom

        elif action_type == "P":
            pass

        else:
            raise "Unknown action type"

        action_plot.scatter(dots)
        action_plot.fig.tight_layout()
        print("Saving large action #{}...".format(action_index))
        action_plot.save(mkpath(output_folder, "history", str(action_index) + ".png"))
        action_plot.close()


if __name__ == "__main__":
    main(query_genome_path, ref_genome_path, sam_file_path, show_plot, output_folder)
