from string import ascii_uppercase
from Bio.SeqIO.FastaIO import SimpleFastaParser
# from matplotlib.patches import ConnectionStyle
from copy import deepcopy
import os

import sys
sys.path.append("src")
from sam_analyze_utils import mkpath, prettifyNumber, distance2, linearApprox, setSettings, Plot


INT_MAX = int(1e9) + 7

# TODO:
# - Left bottom
# - (large7 - left side) - display all dots


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
    "dot_skip_rate": 10,
    "dotsize": 0.1,
    "fontsize": 8,
    "figsize": (10, 7),

    "min_event_size": int(5e3),
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

query_genome_path = "samples/large7/large_genome1.fasta"
ref_genome_path = "samples/large7/large_genome2.fasta"
sam_file_path = "BWA/large7/bwa_output.sam"
show_plot = True
output_folder = "tests/large7"

# query_genome_path = "samples/small/source.fasta"
# ref_genome_path = "samples/small/duplication.fasta"
# sam_file_path = "BWA/small/duplication/bwa_output.sam"
# show_plot = True
# output_folder = "tests/small/duplication"

# \-----TESTING SETTINGS-----/ #


def dotOnLineY(x1, y1, x2, y2, target_x):
    return y1 + (y2 - y1) * ((target_x - x1) / (x2 - x1))


def analyze(query_genome_path: str, ref_genome_path: str, sam_file_path: str, show_plot: bool, output_folder: str, settings: dict):
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
    print("Reference: {} [{}]\n".format(ref_genome_name, ref_genome_length))

    # return
# ====================================================================================================================================================================
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
# ====================================================================================================================================================================
    # Parse CIGAR and create a list of all actions

    m_actions = []
    for position, flags, rid_size, actions in sam_data:

        rotated = lambda ref_pos: ref_genome_length - ref_pos if 4 in flags else ref_pos
        cur_query_pos = position
        cur_ref_pos = 0

        for length, action_type in actions:

            if action_type in ('S', 'H'):
                cur_ref_pos += length

            elif action_type == 'M':
                m_actions.append([cur_query_pos, cur_ref_pos, length, (4 in flags)])
                cur_query_pos += length
                cur_ref_pos += length

            elif action_type == 'I':
                cur_ref_pos += length

            elif action_type == 'D':
                cur_query_pos += length

    # return
# ====================================================================================================================================================================
    # Creating plot
    print("Creating plot...")

    plot = Plot("Main plot", settings["fontsize"], settings["grid_size"], settings["figsize"], query_genome_name, ref_genome_name)
    plot.legendLine({
        "Insertion": "#0f0",
        "Deletion": "#f00",
        "Duplication": "#f0f",
        "Back Duplication": "#0ff"
    }, fontsize=settings["fontsize"], lw=2)

    # return
# ====================================================================================================================================================================
    # Creating dots
    print("Creating dots...", end="")

    graph = [[] for _ in range(query_genome_length + 1)]
    count = 0

    for cur_query_pos, cur_ref_pos, length, rotated in m_actions:
        if rotated:
            cur_ref_pos = ref_genome_length - cur_ref_pos
        for _ in range(length):
            graph[cur_query_pos].append(cur_ref_pos)
            cur_query_pos += 1
            cur_ref_pos += (-1 if rotated else 1)

        count += length

    print(" {}".format(prettifyNumber(count)))  # Can be with optional compress

    # return
# ====================================================================================================================================================================
    # Counting lines
    print("Counting lines...", end="")

    lines_join_size2 = settings["lines_join_size"] ** 2
    line_min_size2 = settings["line_min_size"] ** 2

    lines = []  # Struct: { [start_x, start_y, end_x, end_y, [dots]] }

    for x in range(0, len(graph), settings["dot_skip_rate"]):
        for y in graph[x]:
            for line in lines:

                if distance2(x, y, *line[4][-1]) <= lines_join_size2 and \
                        (len(line[4]) == 1 or distance2(x, y, *line[4][-2]) <= lines_join_size2):
                    line[4].append([x, y])
                    break
            else:
                lines.append([None, None, None, None, [[x, y]]])

    for line in lines:
        line[4].sort()

        line[0], line[1] = line[4][0]
        line[2], line[3] = line[4][-1]

        if len(line[4]) >= 2:
            k, b = linearApprox(line[4])  # \
            line[1] = k * line[0] + b     # |-> Approximation
            line[3] = k * line[2] + b     # /

        # line[4] = line[4][::settings["dot_skip_rate"]]  # Optional compress

    lines = [line for line in lines if distance2(line[0], line[1], line[2], line[3]) >= line_min_size2]

    lines.sort(key=lambda line: (line[0], line[1]))

    # print(*[line[0:4] for line in lines], sep='\n')

    print(" {} lines".format(len(lines)))

    # return
# ====================================================================================================================================================================
    # Rotations
    print("Joining rotations...")

    def rotateLines(lines):
        print("Rotation:")
        rotation_center, rotation_actions = [None] * len(lines), []
        rotation_start = [(i if lines[i][3] < lines[i][1] else None) for i in range(len(lines))]

        last_rotation = None if rotation_start[0] is None else 0
        for i in range(1, len(lines)):
            if rotation_start[i] is not None:
                if rotation_start[i - 1] is not None:
                    rotation_start[i] = rotation_start[i - 1]
                if last_rotation is not None and lines[i][0] <= lines[last_rotation][2] + settings["rotations_join_size"]:
                    rotation_start[i] = rotation_start[last_rotation]

                last_rotation = i

        cur_line = len(lines) - 1
        while cur_line >= 0:
            if rotation_start[cur_line] is None:
                cur_line -= 1
                continue

            min_value, max_value = INT_MAX, -INT_MAX

            for i in range(rotation_start[cur_line], cur_line + 1):
                min_value = min(min_value, lines[i][1], lines[i][3])
                max_value = max(max_value, lines[i][1], lines[i][3])

            for i in range(rotation_start[cur_line], cur_line + 1):
                rotation_center[i] = (min_value + max_value) // 2

            if lines[cur_line][2] - lines[rotation_start[cur_line]][0] >= settings["min_event_size"]:
                rotation_actions.append(["Rotation", rotation_start[cur_line], cur_line, rotation_center[cur_line]])

            # plot.ax.axhline(y=rotation_center[i], color='#ff0', linestyle='-')

            cur_line = rotation_start[cur_line] - 1

        rotation_actions.sort(key=lambda action: (action[1], action[2]))

        print("rotation_start:", rotation_start)
        print("rotation_center:", rotation_center)
        print("rotation_actions:", rotation_actions)

        rotated_lines = deepcopy(lines)  # TODO: remove deepcopy
        for i in range(len(lines)):
            if rotation_center[i] is not None:
                line = rotated_lines[i]
                line[1] -= (line[1] - rotation_center[i]) * 2
                line[3] -= (line[3] - rotation_center[i]) * 2

                for j in range(len(line[4])):
                    line[4][j][1] -= (line[4][j][1] - rotation_center[i]) * 2

        return rotated_lines, rotation_actions

    rotated_lines, rotation_actions = deepcopy(lines), []

    for _ in range(100):
        for i in range(len(rotated_lines)):
            if rotated_lines[i][3] < rotated_lines[i][1]:
                rotated_lines, cur_rotation_actions = rotateLines(rotated_lines)
                # TODO: endless loop?
                rotation_actions += cur_rotation_actions
                break
        else:
            break

    else:
        exit("Endless loop!!!")

    print([line[:4] for line in lines])

    # return
# ====================================================================================================================================================================
    # Handle events (actions)
    print("Handling lines (actions)...")

    large_actions = []

    last_line = rotated_lines[0]
    for line_index in range(1, len(rotated_lines)):
        cur_line = rotated_lines[line_index]
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
                large_actions.append(["Translocation", last_query_end, last_ref_end, back_duplication_length])

            plot.line(cur_query_start, last_ref_end, cur_query_start, cur_ref_start, color="#0ff")

        else:
            print("Unknown!!!")

        if cur_query_end >= last_query_end:
            last_line = cur_line

    large_actions = sorted(large_actions, key=lambda action: -action[3])

    # print(large_actions)

    # return
# ====================================================================================================================================================================
    # Plotting dots and lines
    print("Plotting dots and lines...")

    for line in lines:
        plot.line(line[0], line[1], line[2], line[3], color="#ccc")

    for line in rotated_lines:
        plot.line(line[0], line[1], line[2], line[3], color="#000")

    dots = []
    for x in range(0, len(graph), settings["dot_skip_rate"]):
        dots += ([x, y] for y in graph[x])

    plot.scatter(dots, dotsize=settings["dotsize"], color="#00f")

    print("Saving plot...")
    # plot.tight()
    plot.save(mkpath(output_folder, "sam_analyze.png"))

    if show_plot:
        print("Showing plot...")
        plot.show()

    plot.clear()

    # return
# ====================================================================================================================================================================
    # Make and save history
    print("Making history...", end="")

    if not os.path.exists(mkpath(output_folder, "history")):
        os.mkdir(mkpath(output_folder, "history"))

    for filename in os.listdir(mkpath(output_folder, "history")):
        os.remove(mkpath(output_folder, "history", filename))

    large_actions = [["Pass", 0, 0, 0]] + rotation_actions + large_actions

    print(" {} images\n".format(len(large_actions)))

    print(rotation_actions)
    print(large_actions)

    for action_index in range(len(large_actions)):
        action = large_actions[action_index]
        print(action)
        action_type, start_query_pos, start_ref_pos, length = action[0:4]  # length - направленная длина!!! (уже нет) наверное, но это не точно. Нет!

        if action_type == "Rotation":

            rotation_start, rotation_end, rotation_center = start_query_pos, start_ref_pos, length

            for i in range(rotation_start, rotation_end + 1):
                for j in range(len(lines[i][4])):
                    lines[i][4][j][1] = lines[i][4][j][1] - (lines[i][4][j][1] - rotation_center) * 2

            # for i in range(len(dots)):
            #     if start_query_pos <= dots[i][0] <= start_query_pos + length:
            #         dots[i][1] = dots[i][1] - (dots[i][1] - rotation_center) * 2

        elif action_type == "Deletion":
            for line in rotated_lines:
                for i in range(len(line[4])):
                    if line[4][i][0] >= start_query_pos + length:
                        line[4][i][0] -= length

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos + length:
                    large_actions[i][1] -= length

        elif action_type == "Insertion":

            for line in rotated_lines:
                new_dots = []
                for dot_x, dot_y in line[4]:
                    if dot_x > start_query_pos:
                        dot_y -= length
                    if dot_x != start_query_pos:
                        new_dots.append([dot_x, dot_y])

                line[4] = new_dots

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] -= length

        elif action_type == "Duplication":
            height, line_index = action[4], action[5]

            new_dots = []
            for dot in rotated_lines[line_index][4]:
                if start_query_pos <= dot[0] <= start_query_pos + length:
                    pass
                else:
                    new_dots.append(dot)
            rotated_lines[line_index][4] = new_dots

            for line in rotated_lines:
                for i in range(len(line[4])):
                    if line[4][i][0] >= start_query_pos:
                        line[4][i][1] -= height

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] -= height

        elif action_type == "Translocation":
            for line in rotated_lines:
                new_dots = []
                for dot_x, dot_y in line[4]:
                    if dot_x > start_query_pos:
                        dot_y += length
                    if dot_x != start_query_pos:
                        new_dots.append([dot_x, dot_y])

                line[4] = new_dots

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] += length

        elif action_type == "Pass":
            pass

        else:
            raise "Unknown action type"

        if action_type in ("Pass", "Rotation"):
            # plot.scatter(dots, dotsize=settings["dotsize"], color="#00f")
            for line in lines:
                plot.scatter(line[4], dotsize=settings["dotsize"], color="#00f")
        else:
            # Adjusting axes (bottom):
            bottom = INT_MAX
            for line in rotated_lines:
                for dot_x, dot_y in line[4]:
                    bottom = min(bottom, dot_y)

            for line in rotated_lines:
                for i in range(len(line[4])):
                    line[4][i][1] -= bottom

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] -= bottom

            for line in rotated_lines:
                plot.scatter(line[4], dotsize=settings["dotsize"], color="#00f")

        print("Saving large action #{}...\n".format(action_index))
        plot.tight()
        plot.save(mkpath(
            output_folder,
            "history",
            "{}{}.png".format(
                str(action_index).zfill(3),
                "" if action_type == "Pass" else " (" + action_type + ")"
            )
        ))
        plot.clear()

    del plot


if __name__ == "__main__":
    analyze(query_genome_path, ref_genome_path, sam_file_path, show_plot, output_folder, SETTINGS)
