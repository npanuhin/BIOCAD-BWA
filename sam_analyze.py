from string import ascii_uppercase
from Bio.SeqIO.FastaIO import SimpleFastaParser
from json import load as json_load  # , dump as json_dump
from copy import deepcopy
import os

import sys
sys.path.append("src")
from utils import mkpath, prtNum, distance2, linearApprox, YCoordOnLine, setSettings

from Line import Line
from Plot import Plot

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

with open(mkpath("src", "CIGAR_FLAGS.json"), 'r', encoding="utf-8") as file:
    CIGAR_FLAGS = json_load(file, encoding="utf-8")


# /-----TESTING SETTINGS-----\ #

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

# \-----TESTING SETTINGS-----/ #


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

    print("Query: {} [{}]".format(query_genome_name, prtNum(query_genome_length)))
    print("Reference: {} [{}]\n".format(ref_genome_name, prtNum(ref_genome_length)))

    # return
# ====================================================================================================================================================================
    # Parse CIGAR and create a list of all actions
    print("Reading SAM file...")

    segments = []
    with open(sam_file_path, 'r', encoding="utf-8") as sam_file:
        for line in (line.strip().split() for line in sam_file if not line.strip().startswith("@")):
            # Quality:
            # mapQ = int(line[4])
            # quality = round((10 ** (mapQ / -10)) * 100, 6)

            # Flags:
            flags_bit = int(line[1])
            flags = set()
            for i in range(len(CIGAR_FLAGS) - 1, -1, -1):
                cur_flag, flags_bit = divmod(flags_bit, 2 ** i)
                if cur_flag:
                    flags.add(i)

            # Rid:
            rid_size = len(line[9])
            if rid_size <= settings["min_rid_size"]:
                continue

            # CIGAR:
            actions = []
            buff = ""
            for char in line[5]:
                if char in ascii_uppercase:
                    actions.append([int(buff), char])
                    buff = ""
                else:
                    buff += char

            rotated = lambda ref_pos: ref_genome_length - ref_pos if 4 in flags else ref_pos

            # Start position:
            cur_query_pos = int(line[3])
            cur_ref_pos = 0

            for length, action_type in actions:

                if action_type in ('S', 'H'):
                    cur_ref_pos += length

                elif action_type == 'M':
                    segments.append([cur_query_pos, cur_ref_pos, length, (4 in flags)])
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
        "Translocation": "#0ff"
    }, fontsize=settings["fontsize"], lw=2)

    # return
# ====================================================================================================================================================================
    # Creating dots
    print("Creating dots...", end="")

    graph = [[] for _ in range(query_genome_length + 1)]
    count = 0

    for cur_query_pos, cur_ref_pos, length, rotated in segments:
        if rotated:
            cur_ref_pos = ref_genome_length - cur_ref_pos
        for _ in range(length):
            graph[cur_query_pos].append(cur_ref_pos)
            cur_query_pos += 1
            cur_ref_pos += (-1 if rotated else 1)

        count += length

    del segments

    print(" {}".format(prtNum(count)))  # Can be with optional compress: count // N

    # return
# ====================================================================================================================================================================
    # Counting lines
    print("Counting lines...", end="")

    lines_join_size2 = settings["lines_join_size"] ** 2
    line_min_size2 = settings["line_min_size"] ** 2

    lines = []

    for x in range(0, len(graph), settings["dot_skip_rate"]):
        for y in graph[x]:
            for line in lines:
                if distance2(x, y, *line.dots[-1]) <= lines_join_size2 and \
                        (len(line.dots) == 1 or distance2(x, y, *line.dots[-2]) <= lines_join_size2):
                    line.dots.append([x, y])
                    break
            else:
                lines.append(Line(dots=[[x, y]]))

    for line in lines:
        line.dots.sort()

        line.start_x, line.start_y = line.dots[0]
        line.end_x, line.end_y = line.dots[-1]

        if len(line.dots) >= 2:
            k, b = linearApprox(line.dots)        # \
            line.start_y = k * line.start_x + b   # |--> Approximation
            line.end_y = k * line.end_x + b       # /

        # line[4] = line[4][::settings["dot_skip_rate"]]  # Optional compress

    lines = [line for line in lines if distance2(line.start_x, line.start_y, line.end_x, line.end_y) >= line_min_size2]

    lines.sort(key=lambda line: (line.start_x, line.start_y))

    print(" {} lines".format(len(lines)))
    # print(*[line.coords for line in lines], sep='\n')

    # return
# ====================================================================================================================================================================
    # Rotations
    print("Joining rotations...")

    def rotateLines(lines):
        # print("Rotation:")
        rotation_center, rotation_actions = [None] * len(lines), []
        rotation_start = [(i if lines[i].end_y < lines[i].start_y else None) for i in range(len(lines))]

        # print(*[line.coords for line in lines], sep='\n')

        # Counting rotation_start
        last_rotation = None
        for line_index in range(len(lines)):
            if rotation_start[line_index] is not None:
                if line_index > 0 and rotation_start[line_index - 1] is not None:
                    rotation_start[line_index] = rotation_start[line_index - 1]
                if last_rotation is not None and lines[line_index].start_x <= lines[last_rotation].end_x + settings["rotations_join_size"]:
                    rotation_start[line_index] = rotation_start[last_rotation]

                # if rotation_start[line_index] is not None and min(lines[line_index].start_y, lines[line_index].end_y) >= max(lines[rotation_start[line_index]].start_y, lines[rotation_start[line_index]].end_y):
                #     rotation_start[line_index] = line_index

                last_rotation = line_index

        # Counting rotation_center and rotation_actions
        line_index = len(lines) - 1
        while line_index >= 0:
            if rotation_start[line_index] is None:
                line_index -= 1
                continue

            min_value, max_value = INT_MAX, -INT_MAX

            for i in range(rotation_start[line_index], line_index + 1):
                min_value = min(min_value, lines[i].start_y, lines[i].end_y)
                max_value = max(max_value, lines[i].start_y, lines[i].end_y)

            for i in range(rotation_start[line_index], line_index + 1):
                rotation_center[i] = (min_value + max_value) // 2

            if lines[line_index].end_x - lines[rotation_start[line_index]].start_x >= settings["min_event_size"]:
                rotation_actions.append(["Rotation", rotation_start[line_index], line_index, rotation_center[line_index]])

            # plot.hline(rotation_center[i], color='#ff0', linestyle='-')

            line_index = rotation_start[line_index] - 1

        rotation_actions.sort(key=lambda action: (action[1], action[2]))

        # print("rotation_start:", rotation_start)
        # print("rotation_center:", rotation_center)
        # print("rotation_actions:", rotation_actions)

        for i in range(len(lines)):
            if rotation_center[i] is not None:
                line = lines[i]
                line.start_y -= (line.start_y - rotation_center[i]) * 2
                line.end_y -= (line.end_y - rotation_center[i]) * 2

                for j in range(len(line.dots)):
                    line.dots[j][1] -= (line.dots[j][1] - rotation_center[i]) * 2

        return rotation_actions

    rotated_lines, rotation_actions = deepcopy(lines), []

    for _ in range(100):  # TODO: endless loop => const 100
        for i in range(len(rotated_lines)):
            if rotated_lines[i].end_y < rotated_lines[i].start_y:
                rotation_actions += rotateLines(rotated_lines)
                break
        else:
            break

    else:
        exit("Endless loop!!!")

    print([line.coords for line in lines])

    # return
# ====================================================================================================================================================================
    # Handle events (actions)
    print("Handling lines (actions)...")

    actions = []

    last = rotated_lines[0]
    for line_index in range(1, len(rotated_lines)):
        cur = rotated_lines[line_index]

        if cur.start_x >= last.end_x and cur.start_y >= last.end_y:  # top right
            insertion_length = cur.start_y - last.end_y
            deletion_length = cur.start_x - last.end_x

            actions.append(["Insertion", last.end_x, last.end_y, insertion_length])
            plot.line(last.end_x, last.end_y, last.end_x, cur.start_y, color="#0f0")

            actions.append(["Deletion", last.end_x, last.end_y, deletion_length])
            plot.line(last.end_x, cur.start_y, cur.start_x, cur.start_y, color="#f00")

        elif cur.start_x < last.end_x and cur.start_y >= last.end_y:  # top left
            tmp_dot_y = YCoordOnLine(*last.coords, cur.start_x)
            insertion_length = cur.start_y - cur.end_y
            duplication_length = last.end_x - cur.start_x
            duplication_height = last.end_y - tmp_dot_y

            actions.append(["Insertion", cur.start_x, last.end_y, insertion_length])
            plot.line(cur.start_x, last.end_y, cur.start_x, cur.start_y, color="#0f0")

            actions.append(["Duplication", cur.start_x, tmp_dot_y, duplication_length, duplication_height, line_index - 1])
            plot.poligon([
                (cur.start_x, tmp_dot_y),
                (cur.start_x, last.end_y),
                (last.end_x, last.end_y)
            ], color="#f0f")

        elif cur.start_x >= last.end_x and cur.start_y < last.end_y:  # bottom right
            deletion_length = cur.start_x - last.end_x
            translocation_length = last.end_y - cur.start_y

            actions.append(["Deletion", last.end_x, last.end_y, deletion_length])
            plot.line(last.end_x, last.end_y, cur.start_x, last.end_y, color="#f00")

            actions.append(["Translocation", last.end_x, last.end_y, translocation_length])
            plot.line(cur.start_x, last.end_y, cur.start_x, cur.start_y, color="#0ff")

        else:
            print("\nUnknown!!!\n")

        if cur.end_x >= last.end_x:
            last = cur

    large_actions = [action for action in actions if action[3] >= settings["min_event_size"]]
    large_actions.sort(key=lambda action: -action[3])

    print("\nactions:", actions)
    print("\nlarge_actions:", large_actions)
    # print()

    # return
# ====================================================================================================================================================================
    # Plotting dots and lines
    print("Plotting dots and lines...")

    for line in lines:
        plot.plotLine(line, color="#fa0")

    for line in rotated_lines:
        plot.plotLine(line, color="#000")

    dots = []  # Optional compress
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

    with open(mkpath(output_folder, "history.txt"), 'w', encoding="utf-8") as history_file:
        for action in large_actions:

            if action[0] == "Rotation":
                print("Rotation from {} (Query) to {} (Query)\n".format(
                    prtNum(int(lines[action[1]].start_x)), prtNum(int(lines[action[2]].end_x))
                ), file=history_file)

            elif action[0] == "Insertion":
                print("Insertion of {}-{} (Ref) to {} (Query)\n".format(
                    prtNum(int(action[2])), prtNum(int(action[2] + action[3])), prtNum(int(action[1]))
                ), file=history_file)

            elif action[0] == "Deletion":
                print("Deletion of {}-{} (Query) from {} (Ref)\n".format(
                    prtNum(int(action[1])), prtNum(int(action[1] + action[3])), prtNum(int(action[2]))
                ), file=history_file)

            elif action[0] == "Duplication":
                print("Duplication of {}-{} (Query) {}-{} (Ref)\n".format(
                    prtNum(int(action[1])), prtNum(int(action[1] + action[3])), prtNum(int(action[2])), prtNum(int(action[2] + action[4]))
                ), file=history_file)

            elif action[0] == "Translocation":
                print("Translocation of {}-END (Query) from {} (Ref) to {} (Ref)\n".format(
                    prtNum(int(action[1])), prtNum(int(action[2] - action[3])), prtNum(int(action[2]))
                ), file=history_file)

    print(" {} images\n".format(len(large_actions)))

    print(rotation_actions)
    print(large_actions)

    for action_index, action in enumerate(large_actions):
        action_type, start_query_pos, start_ref_pos, length = action[0:4]

        if action_type == "Rotation":

            rotation_start, rotation_end, rotation_center = start_query_pos, start_ref_pos, length

            for i in range(rotation_start, rotation_end + 1):
                for j in range(len(lines[i].dots)):
                    lines[i].dots[j][1] = lines[i].dots[j][1] - (lines[i].dots[j][1] - rotation_center) * 2

            # for i in range(len(dots)):
            #     if start_query_pos <= dots[i][0] <= start_query_pos + length:
            #         dots[i][1] = dots[i][1] - (dots[i][1] - rotation_center) * 2

        elif action_type == "Insertion":
            for line in rotated_lines:
                new_dots = []
                for dot_x, dot_y in line.dots:
                    if dot_x > start_query_pos:
                        dot_y -= length
                    if dot_x != start_query_pos:
                        new_dots.append([dot_x, dot_y])

                line.dots = new_dots

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] -= length

        elif action_type == "Deletion":
            for line in rotated_lines:
                for i in range(len(line.dots)):
                    if line.dots[i][0] >= start_query_pos + length:
                        line.dots[i][0] -= length

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos + length:
                    large_actions[i][1] -= length

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
                for i in range(len(line.dots)):
                    if line.dots[i][0] >= start_query_pos:
                        line.dots[i][1] -= height

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] -= height

        elif action_type == "Translocation":
            for line in rotated_lines:
                new_dots = []
                for dot_x, dot_y in line.dots:
                    if dot_x > start_query_pos:
                        dot_y += length
                    if dot_x != start_query_pos:
                        new_dots.append([dot_x, dot_y])

                line.dots = new_dots

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] += length

        elif action[0] == "Pass":
            pass

        else:
            raise ValueError("Unknown action type")

        if action[0] in ("Pass", "Rotation"):
            # plot.scatter(dots, dotsize=settings["dotsize"], color="#00f")
            for line in lines:
                plot.scatter(line.dots, dotsize=settings["dotsize"], color="#00f")
        else:
            # Adjusting axes (bottom):
            bottom = INT_MAX
            for line in rotated_lines:
                for dot_x, dot_y in line.dots:
                    bottom = min(bottom, dot_y)

            for line in rotated_lines:
                for i in range(len(line.dots)):
                    line.dots[i][1] -= bottom

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i][1] >= start_query_pos:
                    large_actions[i][2] -= bottom

            for line in rotated_lines:
                plot.scatter(line.dots, dotsize=settings["dotsize"], color="#00f")

        print("Saving large action #{}{}...\n".format(action_index, "" if action[0] == "Pass" else " (" + action[0] + ")"))
        plot.tight()
        plot.save(mkpath(
            output_folder,
            "history",
            "{}{}.png".format(
                str(action_index).zfill(3),
                "" if action[0] == "Pass" else " (" + action[0] + ")"
            )
        ))
        plot.clear()

    del plot


if __name__ == "__main__":
    analyze(query_genome_path, ref_genome_path, sam_file_path, show_plot, output_folder, SETTINGS)
