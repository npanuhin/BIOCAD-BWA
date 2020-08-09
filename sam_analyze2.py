from string import ascii_uppercase
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio.SeqIO.FastaIO import SimpleFastaParser
# from matplotlib.patches import ConnectionStyle
from json import loads as json_loads
import sys
import os

INT_MAX = int(1e9) + 7

# TODO:
# - Fix duplication: insertion must be way smaller


# !!! X - query, Y - ref !!!

# Small
# GRID_SIZE = 1e2
# MIN_RID_SIZE = 1
# DOT_SKIP_RATE = 1
# DOT_SIZE = 0.1
# MIN_EVENT_SIZE = 3
# ROTATION_JOIN_SIZE = 1e1
# LINES_JOIN_SIZE = 1e1
# LINE_MIN_SIZE = 1e1

# Large
GRID_SIZE = int(1e5)
MIN_RID_SIZE = int(1e3)
DOT_SKIP_RATE = 10
DOT_SIZE = 0.1
MIN_EVENT_SIZE = int(1e3)
ROTATION_JOIN_SIZE = int(1e5)
LINES_JOIN_SIZE = int(1e3)
LINE_MIN_SIZE = int(1e1)


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

# query_genome_path = "samples/small/source.fasta"
# ref_genome_path = "samples/small/duplication.fasta"
# sam_file_path = "BWA/small/duplication/bwa_output.sam"
# show_plot = True
# output_folder = "tests/small/duplication"

# ---SETTINGS--- #


def mkpath(*paths):
    return os.path.normpath(os.path.join(*paths))


def prettifyNumber(num):
    return "{:,}".format(num)


# def equalE(value1, value2, E):
#     '''Epsilon comparison'''
#     return value1 - E < value2 < value1 + E


def distance2(x1, y1, x2, y2):
    return (x1 - x2) ** 2 + (y1 - y2) ** 2


def dotOnLineY(x1, y1, x2, y2, target_x):
    return y1 + (y2 - y1) * ((target_x - x1) / (x2 - x1))


def setSettings(settings_path):
    global GRID_SIZE, MIN_RID_SIZE, DOT_SKIP_RATE, DOT_SIZE, MIN_EVENT_SIZE, ROTATION_JOIN_SIZE, \
        LINES_JOIN_SIZE, LINE_MIN_SIZE, FIGSIZE, FONT_SIZE, CIGAR_FLAGS

    with open(settings_path, 'r', encoding="utf-8") as settings_file:
        settings = json_loads(settings_file.read().strip())

    for key in settings:
        if isinstance(getattr(sys.modules[__name__], key), int):
            setattr(sys.modules[__name__], key, int(settings[key]))

        elif isinstance(getattr(sys.modules[__name__], key), float):
            setattr(sys.modules[__name__], key, float(settings[key]))

        else:
            setattr(sys.modules[__name__], key, settings[key])


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

    def legendLine(self, legend_dict, fontsize=FONT_SIZE, *line_args, **line_kwargs):
        legend_objects, legend_names = [], []
        for name in legend_dict:
            legend_names.append(name)
            legend_objects.append(Line2D([0], [0], color=legend_dict[name], *line_args, **line_kwargs))

        self.ax.legend(legend_objects, legend_names, fontsize=fontsize)

    def save(self, path, *args, **kwargs):
        self.fig.savefig(path, dpi=400, *args, **kwargs)

    def close(self):
        plt.close(self.fig)


def main(query_genome_path, ref_genome_path, sam_file_path, show_plot, output_folder):

    if os.path.exists(mkpath(output_folder, "settings.json")):
        setSettings(mkpath(output_folder, "settings.json"))

    with open(query_genome_path, 'r', encoding="utf-8") as file:
        for name, sequence in SimpleFastaParser(file):
            query_genome_name = name
            query_genome_length = len(sequence)
            break

    with open(ref_genome_path, 'r', encoding="utf-8") as file:
        for name, sequence in SimpleFastaParser(file):
            ref_genome_name = name
            ref_genome_length = len(sequence)
            break

    print("Query: {} [{}]".format(query_genome_name, query_genome_length))
    print("Reference: {} [{}]".format(ref_genome_name, ref_genome_length))
    print()

    sam_data = []

    with open(sam_file_path, 'r', encoding="utf-8") as sam_file:

        for line in (line.strip().split() for line in sam_file if not line.strip().startswith("@")):
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
            rid_size = len(line[9])
            if rid_size <= MIN_RID_SIZE:
                continue

            # CIGAR
            cigar = line[5]
            actions = []
            buff = ""
            for char in cigar:
                if char in ascii_uppercase:
                    actions.append([int(buff), char])
                    buff = ""
                else:
                    buff += char

            sam_data.append([position, flags, rid_size, actions])

    # return
# ========================================================================================================================================

    bwa_actions = []

    for position, flags, rid_size, actions in sam_data:

        start_query_pos, start_ref_pos, end_query_pos, end_ref_pos = INT_MAX, INT_MAX, 0, 0
        cur_query_pos = position
        cur_ref_pos = 0

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

        print("{} -> {}".format(prettifyNumber(start_query_pos), prettifyNumber(end_query_pos)))
        # for flag in flags:
        #     if flag == 11:  # Disable flag №11
        #         continue
        #     print(CIGAR_FLAGS[flag])

        # for length, action_type in actions:
        #     print("{}-{}|".format(length, action_type), end="")
        # print()

        rotation_center = ((start_ref_pos + end_ref_pos) / 2) if 4 in flags else None

        cur_query_pos = position
        cur_ref_pos = 0

        for length, action_type in actions:

            if action_type not in ('S', 'H'):
                bwa_actions.append([cur_query_pos, cur_ref_pos, length, action_type, rotation_center])

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

    # return
# ========================================================================================================================================
    # Join rotations for bwa_actions + create rotations for large_actions
    # TODO: not count I and D actions in ref/query end

    rotations = []

    def getActionQueryEnd(action):
        if action[3] in ('M', 'I'):
            return action[0] + action[2]
        return action[0]

    def getActionRefEnd(action):
        if action[3] in ('M', 'D'):
            return action[1] + action[2]
        return action[1]

    bwa_actions.sort(key=lambda action: action[0])

    rotation_block_start, rotation_block_end = 0, 0
    last_rotation_query_end = None
    for i in range(len(bwa_actions)):
        action = bwa_actions[i]
        start_query_pos, start_ref_pos, length, action_type, rotation_center = action
        if rotation_center is None:
            continue

        if last_rotation_query_end is None:
            rotation_block_start = i

        elif last_rotation_query_end + ROTATION_JOIN_SIZE >= start_query_pos:
            pass

        else:
            new_rotation_center = (bwa_actions[rotation_block_start][1] + getActionRefEnd(bwa_actions[rotation_block_end])) / 2
            for j in range(rotation_block_start, rotation_block_end + 1):
                bwa_actions[j][4] = new_rotation_center

            block_length = getActionQueryEnd(bwa_actions[rotation_block_end]) - bwa_actions[rotation_block_start][0]
            if block_length >= MIN_EVENT_SIZE:
                rotations.append(["Rotation", bwa_actions[rotation_block_start][0], bwa_actions[rotation_block_start][1], block_length, ref_genome_length - new_rotation_center])

            rotation_block_start = i

        rotation_block_end = i
        last_rotation_query_end = getActionQueryEnd(action)

    if last_rotation_query_end is not None:
        new_rotation_center = (bwa_actions[rotation_block_start][1] + getActionRefEnd(bwa_actions[rotation_block_end])) / 2
        for j in range(rotation_block_start, rotation_block_end + 1):
            bwa_actions[j][4] = new_rotation_center

        block_length = getActionQueryEnd(bwa_actions[rotation_block_end]) - bwa_actions[rotation_block_start][0]
        if block_length >= MIN_EVENT_SIZE:
            rotations.append(["Rotation", bwa_actions[rotation_block_start][0], bwa_actions[rotation_block_start][1], block_length, ref_genome_length - new_rotation_center])

    # return
# ========================================================================================================================================
    # Create dots

    bwa_actions.sort(key=lambda action: action[0])

    plot = Plot("Main", query_genome_name, ref_genome_name)
    plot.legendLine({
        "Insertion": "#0f0",
        "Deletion": "#f00",
        "Duplication": "#f0f",
        "Back Duplication": "#0ff"
    }, lw=2)

    last_query_end, last_ref_end = None, None
    dots, ghost_dots, rotated_dots = [], [], []

    graph = [[] for _ in range(query_genome_length + 1)]

    for action_index in range(len(bwa_actions)):

        start_query_pos, start_ref_pos, length, action_type, rotation_center = bwa_actions[action_index]
        cur_query_pos, cur_ref_pos = start_query_pos, start_ref_pos

        if rotation_center is None:
            rotated = lambda cur_ref_pos: cur_ref_pos
        else:
            rotated = lambda cur_ref_pos: (ref_genome_length - cur_ref_pos) + ((ref_genome_length - rotation_center) - (ref_genome_length - cur_ref_pos)) * 2

        if action_type == 'M':
            for i in range(length):
                if rotation_center is None:
                    dots.append([cur_query_pos, cur_ref_pos])
                else:
                    dots.append([cur_query_pos, ref_genome_length - cur_ref_pos])
                    ghost_dots.append([cur_query_pos, rotated(cur_ref_pos)])

                rotated_dots.append([cur_query_pos, rotated(cur_ref_pos)])
                graph[cur_query_pos].append(int(rotated(cur_ref_pos)))

                cur_query_pos += 1
                cur_ref_pos += 1

        elif action_type == 'I':
            cur_ref_pos += length

        elif action_type == 'D':
            cur_query_pos += length

        if last_query_end is None or cur_query_pos >= last_query_end:
            last_query_end = cur_query_pos

    dots = dots[::DOT_SKIP_RATE]
    ghost_dots = ghost_dots[::DOT_SKIP_RATE]
    rotated_dots = rotated_dots[::DOT_SKIP_RATE]

    print("Dots count: {} + {}".format(prettifyNumber(len(dots)), prettifyNumber(len(ghost_dots))))

    # return
# ========================================================================================================================================
    # Count lines

    print("Counting lines...")

    lines = []

    LINES_JOIN_SIZE2 = LINES_JOIN_SIZE ** 2
    LINE_MIN_SIZE2 = LINE_MIN_SIZE ** 2

    for dot_x in range(0, len(graph), DOT_SKIP_RATE):
        for dot_y in graph[dot_x]:
            for i in range(len(lines)):

                if distance2(dot_x, dot_y, lines[i][2], lines[i][3]) <= LINES_JOIN_SIZE2:
                    lines[i][0] = min(lines[i][0], dot_x)
                    lines[i][1] = min(lines[i][1], dot_y)
                    lines[i][2] = max(lines[i][2], dot_x)
                    lines[i][3] = max(lines[i][3], dot_y)
                    lines[i][4].append([dot_x, dot_y])
                    break
            else:
                lines.append([dot_x, dot_y, dot_x, dot_y, [[dot_x, dot_y]]])

    lines = [line for line in lines if distance2(line[0], line[1], line[2], line[3]) >= LINE_MIN_SIZE2]

    print("{} lines".format(len(lines)))

    # return
# ========================================================================================================================================
    # Handle events

    large_actions = []

    lines.sort(key=lambda line: line[0])

    for line in lines:
        plot.line(line[0], line[1], line[2], line[3], color="black")

    last_line = lines[0]
    for line_index in range(1, len(lines)):
        cur_line = lines[line_index]
        last_query_start, last_ref_start, last_query_end, last_ref_end, last_dots = last_line
        cur_query_start, cur_ref_start, cur_query_end, cur_ref_end, cur_dots = cur_line

        if cur_query_start >= last_query_end and cur_ref_start >= last_ref_end:  # top right
            insertion_length = cur_ref_start - last_ref_end
            deletion_length = cur_query_start - last_query_end

            if insertion_length >= MIN_EVENT_SIZE:
                large_actions.append(["Insertion", last_query_end, last_ref_end, insertion_length])

            plot.line(last_query_end, last_ref_end, last_query_end, cur_ref_start, color="#0f0")

            if deletion_length >= MIN_EVENT_SIZE:
                large_actions.append(["Deletion", last_query_end, last_ref_end, deletion_length])

            plot.line(last_query_end, cur_ref_start, cur_query_start, cur_ref_start, color="#f00")

        elif cur_query_start < last_query_end and cur_ref_start >= last_ref_end:  # top left
            insertion_length = dotOnLineY(*cur_line[:4], last_query_end) - last_ref_end
            duplication_length = last_query_end - cur_query_start

            if insertion_length >= MIN_EVENT_SIZE:
                large_actions.append(["Insertion", cur_query_start, cur_ref_start - insertion_length, insertion_length])

            plot.line(cur_query_start, cur_ref_start - insertion_length, cur_query_start, cur_ref_start, color="#0f0")

            if duplication_length >= MIN_EVENT_SIZE:
                large_actions.append(["Duplication", cur_query_start, cur_ref_start, duplication_length, line_index - 1])

            plot.line(cur_query_start, cur_ref_start - insertion_length, last_query_end, cur_ref_start - insertion_length, color="#f0f")

        elif cur_query_start >= last_query_end and cur_ref_start < last_ref_end:  # bottom right
            deletion_length = cur_query_start - last_query_end
            back_duplication_length = last_ref_end - cur_ref_start

            if deletion_length >= MIN_EVENT_SIZE:
                large_actions.append(["Deletion", last_query_end, last_ref_end, deletion_length])

            plot.line(last_query_end, last_ref_end, cur_query_start, last_ref_end, color="#f00")

            if back_duplication_length >= MIN_EVENT_SIZE:
                large_actions.append(["Back Dupication", last_query_end, last_ref_end, back_duplication_length])

            plot.line(cur_query_start, last_ref_end, cur_query_start, cur_ref_start, color="#0ff")

        else:
            print("Unknown!!!")

        if cur_query_end >= last_query_end:
            last_line = cur_line

    large_actions = rotations + sorted(large_actions, key=lambda action: -action[3])

    print(large_actions)

    # return
# ========================================================================================================================================
    # Save and show main plot

    plot.fig.tight_layout()
    print("Saving plot...\n")
    plot.save(mkpath(output_folder, "sam_analyze.png"))

    if show_plot:
        print("Showing plot...\n")
        plt.show()

    del plot

    plot = Plot("Main (dots)", query_genome_name, ref_genome_name)

    plot.scatter(dots)
    if ghost_dots:
        plot.scatter(ghost_dots, color="#ccc")

    plot.fig.tight_layout()
    print("Saving dot plot...\n")
    plot.save(mkpath(output_folder, "sam_analyze (dot plot).png"))

    del plot

    # return
# ========================================================================================================================================
    # Save history

    if not os.path.exists(mkpath(output_folder, "history")):
        os.mkdir(mkpath(output_folder, "history"))

    for filename in os.listdir(mkpath(output_folder, "history")):
        os.remove(mkpath(output_folder, "history", filename))

    large_actions = [["Pass", 0, 0, 0]] + large_actions

    print("History size: {} images\n".format(len(large_actions)))

    for action_index in range(len(large_actions)):
        action = large_actions[action_index]
        print(action)
        action_type, start_query_pos, start_ref_pos, length = action[0:4]  # length - направленная длина!!! (уже нет) наверное, но это не точно

        action_plot = Plot("large_action{}".format(action_index + 1), nameX=query_genome_name, nameY=ref_genome_name)

        if action_type == "Rotation":

            rotation_center = action[4]

            # for line in lines:
            #     for i in range(len(line[4])):
            #         if start_query_pos <= line[4][i][0] <= start_query_pos + length:
            #             line[4][i][1] = line[4][i][1] - (line[4][i][1] - rotation_center) * 2

            for i in range(len(dots)):
                if start_query_pos <= dots[i][0] <= start_query_pos + length:
                    dots[i][1] = dots[i][1] - (dots[i][1] - rotation_center) * 2

        elif action_type == "Deletion":

            for line in lines:
                for i in range(len(line[4])):
                    if line[4][i][0] >= start_query_pos + length:
                        line[4][i][0] -= length

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos + length:
                    large_actions[i][1] -= length

        elif action_type == "Insertion":

            for line in lines:
                for i in range(len(line[4])):
                    if line[4][i][0] >= start_query_pos:
                        line[4][i][1] -= length

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] -= length

        elif action_type == "Duplication":
            line_index = action[4]

            new_dots = []

            for dot in lines[line_index][4]:
                if start_query_pos <= dot[0] <= start_query_pos + length:
                    pass
                else:
                    new_dots.append(dot)

            lines[line_index][4] = new_dots

        elif action_type == "Back Dupication":
            pass

        elif action_type == "Pass":
            pass

        else:
            raise "Unknown action type"

        if action_type in ("Pass", "Rotation"):
            action_plot.scatter(dots, color="blue")
        else:

            # Adjusting axes (bottom):
            bottom = INT_MAX
            for line in lines:
                for dot_x, dot_y in line[4]:
                    bottom = min(bottom, dot_y)

            for line in lines:
                for i in range(len(line[4])):
                    line[4][i][1] -= bottom

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] -= bottom

            # Plotting:

            for line in lines:
                action_plot.scatter(line[4], color="blue")

        action_plot.fig.tight_layout()
        print("Saving large action #{}...\n".format(action_index))
        action_plot.save(mkpath(output_folder, "history", str(action_index).zfill(3) + ".png"))
        action_plot.close()


if __name__ == "__main__":
    main(query_genome_path, ref_genome_path, sam_file_path, show_plot, output_folder)
