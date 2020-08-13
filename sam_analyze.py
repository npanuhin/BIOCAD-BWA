from string import ascii_uppercase
from Bio.SeqIO.FastaIO import SimpleFastaParser
# from matplotlib.patches import ConnectionStyle
from json import load as json_load
import os

import sys
sys.path.append("src")
from maplotlib_wrapper import Plot


INT_MAX = int(1e9) + 7

# TODO:
# - Not count I and D actions in ref/query end in main funtion (rotation)
# - Left bottom


# !!! X - query, Y - ref !!!

# Small
# SETTINGS = {
#     "grid_size": 100,
#     "min_rid_size": 1,
#     "dot_skip_rate": 1,
#     "dotsize": 0.1,
#     "fontsize": 10,
#     "figsize": (10, 7),

#     "min_event_size": 3,
#     "rotations_join_size": 10,
#     "lines_join_size": 5,
#     "line_min_size": 10
# }

# Large
SETTINGS = {
    "grid_size": int(1e5),
    "min_rid_size": int(1e3),
    "dot_skip_rate": 500,
    "dotsize": 0.1,
    "fontsize": 10,
    "figsize": (10, 7),

    "min_event_size": int(1e3),
    "rotations_join_size": int(1e5),
    "lines_join_size": "$min_event_size + 3",
    "line_min_size": "$min_event_size"
}

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


# /-----TESTING SETTINGS-----\ #

query_genome_path = "samples/large3/large_genome1.fasta"
ref_genome_path = "samples/large3/large_genome2.fasta"
sam_file_path = "BWA/large3/bwa_output.sam"
show_plot = False
output_folder = "tests/large3"

# query_genome_path = "samples/small/source.fasta"
# ref_genome_path = "samples/small/deletion.fasta"
# sam_file_path = "BWA/small/deletion/bwa_output.sam"
# show_plot = True
# output_folder = "tests/small/deletion"

# \-----TESTING SETTINGS-----/ #


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


def setSettings(settings, alternative_settings_path=None):
    if alternative_settings_path is not None and os.path.exists(alternative_settings_path):
        with open(alternative_settings_path, 'r', encoding="utf-8") as settings_file:
            settings.update(json_load(settings_file))

    counted = {key: settings[key] for key in settings if not isinstance(settings[key], str) or settings[key][0] != '$'}

    while True:
        for key in settings:
            if key not in counted:
                try:
                    settings[key] = eval(settings[key][1:], None, counted)
                except NameError:
                    continue
                counted[key] = settings[key]
                break
        else:
            break


def linearApprox(dots):
    n, sumx, sumy, sumx2, sumxy = len(dots), 0, 0, 0, 0
    for x, y in dots:
        sumx += x
        sumy += y
        sumx2 += x ** 2
        sumxy += x * y

    k = (n * sumxy - (sumx * sumy)) / (n * sumx2 - sumx * sumx)
    b = (sumy - k * sumx) / n
    return k, b


def main(query_genome_path: str, ref_genome_path: str, sam_file_path: str, show_plot: bool, output_folder: str, settings: dict):

    print("---| {} |---".format(output_folder))

    setSettings(settings, mkpath(output_folder, "settings.json"))

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

    # return
# ========================================================================================================================================
    # Read SAM file
    print("Reading SAM file...")

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
            if rid_size <= settings["min_rid_size"]:
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
    # Parse CIGAR and create a list of all actions

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

        # print("{} -> {}".format(prettifyNumber(start_query_pos), prettifyNumber(end_query_pos)))
        # for flag in flags:  # Print flags
        #     if flag == 11:  # Disable flag №11
        #         continue
        #     print(CIGAR_FLAGS[flag])

        # for length, action_type in actions:  # Print CIGAR
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
    print("Joining rotations...")

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

        elif last_rotation_query_end + settings["rotations_join_size"] >= start_query_pos:
            pass

        else:
            new_rotation_center = (bwa_actions[rotation_block_start][1] + getActionRefEnd(bwa_actions[rotation_block_end])) / 2
            for j in range(rotation_block_start, rotation_block_end + 1):
                bwa_actions[j][4] = new_rotation_center

            block_length = getActionQueryEnd(bwa_actions[rotation_block_end]) - bwa_actions[rotation_block_start][0]
            if block_length >= settings["min_event_size"]:
                rotations.append(["Rotation", bwa_actions[rotation_block_start][0], bwa_actions[rotation_block_start][1], block_length, ref_genome_length - new_rotation_center])

            rotation_block_start = i

        rotation_block_end = i
        last_rotation_query_end = getActionQueryEnd(action)

    if last_rotation_query_end is not None:
        new_rotation_center = (bwa_actions[rotation_block_start][1] + getActionRefEnd(bwa_actions[rotation_block_end])) / 2
        for j in range(rotation_block_start, rotation_block_end + 1):
            bwa_actions[j][4] = new_rotation_center

        block_length = getActionQueryEnd(bwa_actions[rotation_block_end]) - bwa_actions[rotation_block_start][0]
        if block_length >= settings["min_event_size"]:
            rotations.append(["Rotation", bwa_actions[rotation_block_start][0], bwa_actions[rotation_block_start][1], block_length, ref_genome_length - new_rotation_center])

    bwa_actions.sort(key=lambda action: action[0])

    # return
# ========================================================================================================================================
    # Create dots
    print("Creating gots...", end="")

    plot = Plot("Main", settings["fontsize"], settings["grid_size"], settings["figsize"], query_genome_name, ref_genome_name)
    plot.legendLine({
        "Insertion": "#0f0",
        "Deletion": "#f00",
        "Duplication": "#f0f",
        "Back Duplication": "#0ff"
    }, fontsize=settings["fontsize"], lw=2)

    last_query_end, last_ref_end = None, None
    dots, ghost_dots = [], []
    graph = [[] for _ in range(query_genome_length + 1)]

    for action_index in range(len(bwa_actions)):
        cur_query_pos, cur_ref_pos, length, action_type, rotation_center = bwa_actions[action_index]

        if rotation_center is None:
            rotated = lambda cur_ref_pos: cur_ref_pos
        else:
            rotated = lambda cur_ref_pos: (ref_genome_length - cur_ref_pos) + ((ref_genome_length - rotation_center) - (ref_genome_length - cur_ref_pos)) * 2

        if action_type == 'M':
            for _ in range(length):
                if rotation_center is None:
                    dots.append([cur_query_pos, cur_ref_pos])
                else:
                    dots.append([cur_query_pos, ref_genome_length - cur_ref_pos])
                    ghost_dots.append([cur_query_pos, rotated(cur_ref_pos)])

                graph[cur_query_pos].append(int(rotated(cur_ref_pos)))

                cur_query_pos += 1
                cur_ref_pos += 1

        elif action_type == 'I':
            cur_ref_pos += length

        elif action_type == 'D':
            cur_query_pos += length

        if last_query_end is None or cur_query_pos >= last_query_end:
            last_query_end = cur_query_pos

    print(" {} + {}".format(prettifyNumber(len(dots)), prettifyNumber(len(ghost_dots))))

    # return
# ========================================================================================================================================
    # Count lines
    print("Counting lines...", end="")

    lines_join_size2 = settings["lines_join_size"] ** 2
    line_min_size2 = settings["line_min_size"] ** 2

    lines = []  # Struct: { [start_x, start_y, end_x, end_y, [dots]] }

    for x in range(0, len(graph)):
        for y in graph[x]:
            for line in lines:

                if distance2(x, y, *line[4][-1]) <= lines_join_size2 and \
                        (len(line[4]) == 1 or distance2(x, y, *line[4][-2]) <= lines_join_size2):
                    line[4].append([x, y])
                    break
            else:
                lines.append([None, None, None, None, [[x, y]]])

    for line in lines:
        line[4].sort(key=lambda dot: dot[0])
        k, b = linearApprox(line[4])

        line[0] = line[4][0][0]
        line[1] = k * line[0] + b

        line[2] = line[4][-1][0]
        line[3] = k * line[2] + b

    lines = [line for line in lines if distance2(line[0], line[1], line[2], line[3]) >= line_min_size2]

    print(" {} lines".format(len(lines)))

    # return
# ========================================================================================================================================
    # Handle events (actions)
    print("Handling lines (actions)...")

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

            if insertion_length >= settings["min_event_size"]:
                large_actions.append(["Insertion", last_query_end, last_ref_end, insertion_length])

            plot.line(last_query_end, last_ref_end, last_query_end, cur_ref_start, color="#0f0")

            if deletion_length >= settings["min_event_size"]:
                large_actions.append(["Deletion", last_query_end, last_ref_end, deletion_length])

            plot.line(last_query_end, cur_ref_start, cur_query_start, cur_ref_start, color="#f00")

        elif cur_query_start < last_query_end and cur_ref_start >= last_ref_end:  # top left
            tmp_dot_y = dotOnLineY(*last_line[:4], cur_query_start)
            insertion_length = cur_ref_start - last_ref_end
            duplication_length = last_query_end - cur_query_start
            duplication_height = last_ref_end - tmp_dot_y

            if insertion_length >= settings["min_event_size"]:
                large_actions.append(["Insertion", cur_query_start, last_ref_end, insertion_length])

            plot.line(cur_query_start, last_ref_end, cur_query_start, cur_ref_start, color="#0f0")

            if duplication_length >= settings["min_event_size"]:
                large_actions.append(["Duplication", cur_query_start, tmp_dot_y, duplication_length, duplication_height, line_index - 1])

            plot.poligon([
                (cur_query_start, tmp_dot_y),
                (cur_query_start, last_ref_end),
                (last_query_end, last_ref_end)
            ], color="#f0f")

        elif cur_query_start >= last_query_end and cur_ref_start < last_ref_end:  # bottom right
            deletion_length = cur_query_start - last_query_end
            back_duplication_length = last_ref_end - cur_ref_start

            if deletion_length >= settings["min_event_size"]:
                large_actions.append(["Deletion", last_query_end, last_ref_end, deletion_length])

            plot.line(last_query_end, last_ref_end, cur_query_start, last_ref_end, color="#f00")

            if back_duplication_length >= settings["min_event_size"]:
                large_actions.append(["Back Dupication", last_query_end, last_ref_end, back_duplication_length])

            plot.line(cur_query_start, last_ref_end, cur_query_start, cur_ref_start, color="#0ff")

        else:
            print("Unknown!!!")

        if cur_query_end >= last_query_end:
            last_line = cur_line

    large_actions = rotations + sorted(large_actions, key=lambda action: -action[3])

    # print(large_actions)

    # return
# ========================================================================================================================================
    print("Compressing dots...")

    dots = dots[::settings["dot_skip_rate"]]
    ghost_dots = ghost_dots[::settings["dot_skip_rate"]]
    for line in lines:
        line[4] = line[4][::settings["dot_skip_rate"]]

    # return
# ========================================================================================================================================
    # Save and show main plot

    print("Saving plot...")
    plot.tight()
    plot.save(mkpath(output_folder, "sam_analyze.png"))

    if show_plot:
        print("Showing plot...")
        plot.show()

    del plot

    plot = Plot("Main (dots)", settings["fontsize"], settings["grid_size"], settings["figsize"], query_genome_name, ref_genome_name)
    plot.scatter(dots, dotsize=settings["dotsize"])
    if ghost_dots:
        plot.scatter(ghost_dots, dotsize=settings["dotsize"], color="#ccc")

    print("Saving dot plot...")
    plot.tight()
    plot.save(mkpath(output_folder, "sam_analyze (dot plot).png"))
    del plot

    # return
# ========================================================================================================================================
    # Save history
    print("Making history...", end="")

    if not os.path.exists(mkpath(output_folder, "history")):
        os.mkdir(mkpath(output_folder, "history"))

    for filename in os.listdir(mkpath(output_folder, "history")):
        os.remove(mkpath(output_folder, "history", filename))

    large_actions = [["Pass", 0, 0, 0]] + large_actions

    print(" {} images\n".format(len(large_actions)))

    history_plot = Plot("history", settings["fontsize"], settings["grid_size"], settings["figsize"], query_genome_name, ref_genome_name)

    for action_index in range(len(large_actions)):
        action = large_actions[action_index]
        print(action)
        action_type, start_query_pos, start_ref_pos, length = action[0:4]  # length - направленная длина!!! (уже нет) наверное, но это не точно

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
            height, line_index = action[4], action[5]

            new_dots = []
            for dot in lines[line_index][4]:
                if start_query_pos <= dot[0] <= start_query_pos + length:
                    pass
                else:
                    new_dots.append(dot)
            lines[line_index][4] = new_dots

            for line in lines:
                for i in range(len(line[4])):
                    if line[4][i][0] >= start_query_pos:
                        line[4][i][1] -= height

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] -= height

        elif action_type == "Back Dupication":
            pass

        elif action_type == "Pass":
            pass

        else:
            raise "Unknown action type"

        if action_type in ("Pass", "Rotation"):
            history_plot.scatter(dots, dotsize=settings["dotsize"], color="blue")
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

            for line in lines:
                history_plot.scatter(line[4], dotsize=settings["dotsize"], color="blue")

        print("Saving large action #{}...\n".format(action_index))
        history_plot.tight()
        history_plot.save(mkpath(output_folder, "history", str(action_index).zfill(3) + ".png"))
        history_plot.clear()


if __name__ == "__main__":
    main(query_genome_path, ref_genome_path, sam_file_path, show_plot, output_folder, SETTINGS)
