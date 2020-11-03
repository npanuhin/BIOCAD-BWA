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
from events import Rotation, Insertion, Deletion, Translocation, Duplication, Pass

# INT_MAX = int(1e9) + 7

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
    "lines_join_size": "$min_event_size + 3",
    "line_min_size": "$min_event_size"
}

with open(mkpath("src", "CIGAR_FLAGS.json"), 'r', encoding="utf-8") as file:
    CIGAR_FLAGS = json_load(file, encoding="utf-8")


# /-----TESTING SETTINGS-----\ #

query_genome_path = "samples/large2/large_genome1.fasta"
ref_genome_path = "samples/large2/large_genome2.fasta"
sam_file_path = "BWA/large2/bwa_output.sam"
show_plot = True
output_folder = "tests/large2"

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
            k, b = linearApprox(line.dots)             # \
            line.start_y = int(k * line.start_x + b)   # |--> Approximation  TODO: int
            line.end_y = int(k * line.end_x + b)       # /

        # line[4] = line[4][::settings["dot_skip_rate"]]  # Optional compress

    lines = [line for line in lines if distance2(line.start_x, line.start_y, line.end_x, line.end_y) >= line_min_size2]

    lines.sort(key=lambda line: (line.start_x, line.start_y))

    print(" {} lines".format(len(lines)))
    print(*lines, sep='\n')

    # return
# ====================================================================================================================================================================
    # Rotations
    print("Counting rotations...")

    def countMetric(lines):
        result = 0
        for line in lines:
            result += abs(line.start_x - line.start_y) ** 2 + abs(line.end_x - line.end_y) ** 2
        return result

    def countMetricWithRotation(lines, rotation, apply_changes=False):
        # print(*([line.start_y, line.end_y] for line in lines))

        rotation_center = (
            min(
                lines[rotation.start_line].start_y, lines[rotation.start_line].end_y,
                lines[rotation.end_line].start_y, lines[rotation.end_line].end_y
            ) + max(
                lines[rotation.start_line].start_y, lines[rotation.start_line].end_y,
                lines[rotation.end_line].start_y, lines[rotation.end_line].end_y
            )
        ) // 2

        for line_index in range(rotation.start_line, rotation.end_line + 1):
            line = lines[line_index]

            line.start_y -= (line.start_y - rotation_center) * 2
            line.end_y -= (line.end_y - rotation_center) * 2

        result = countMetric(lines)

        if not apply_changes:
            for line_index in range(rotation.start_line, rotation.end_line + 1):
                line = lines[line_index]

                line.start_y -= (line.start_y - rotation_center) * 2
                line.end_y -= (line.end_y - rotation_center) * 2

        if apply_changes:
            rotation.rotation_center = rotation_center
            for line_index in range(rotation.start_line, rotation.end_line + 1):
                line = lines[line_index]
                for j in range(len(line.dots)):
                    line.dots[j][1] -= (line.dots[j][1] - rotation_center) * 2

        return result

    possible_rotations = []

    for start_line in range(len(lines)):
        for end_line in range(start_line, len(lines)):
            rotation_center = (
                min(
                    lines[start_line].start_y, lines[start_line].end_y,
                    lines[end_line].start_y, lines[end_line].end_y
                ) + max(
                    lines[start_line].start_y, lines[start_line].end_y,
                    lines[end_line].start_y, lines[end_line].end_y
                )
            ) // 2
            possible_rotations.append(Rotation(start_line, end_line, rotation_center))

    cur_metric_value = countMetric(lines)

    rotated_lines = deepcopy(lines)
    rotation_actions = []

    # print("\nPossible rotations:", *possible_rotations, sep='\n')

    while True:
        best_metric_value = float('inf')
        best_rotation_index = 0

        for i, rotation in enumerate(possible_rotations):

            min_line_center, max_line_center = float('inf'), float("-inf")
            for line_index in range(rotation.start_line, rotation.end_line + 1):
                min_line_center = min(min_line_center, rotated_lines[line_index].center_y)
                max_line_center = max(max_line_center, rotated_lines[line_index].center_y)


            if cur_metric < best_metric_value:
                best_metric_value = cur_metric
                best_rotation_index = i

        if best_metric_value >= cur_metric_value:
            break

        # print("\n{} -> {}".format(possible_rotations[best_rotation_index], cur_metric_value))

        cur_metric_value = countMetricWithRotation(rotated_lines, possible_rotations[best_rotation_index], apply_changes=True)

        rotation_actions.append(possible_rotations[best_rotation_index])

    print("\nRotation actions:", *rotation_actions, sep='\n')
    print("\nRotated lines:", *rotated_lines, sep='\n')

    # return
# ====================================================================================================================================================================
    # Handle events
    print("Handling events...")

    actions = []

    last = rotated_lines[0]
    for line_index in range(1, len(rotated_lines)):
        cur = rotated_lines[line_index]

        if cur.start_x >= last.end_x and cur.start_y >= last.end_y:  # top right
            insertion_length = cur.start_y - last.end_y
            deletion_length = cur.start_x - last.end_x

            actions.append(Insertion(last.end_x, last.end_y, insertion_length))
            plot.line(last.end_x, last.end_y, last.end_x, cur.start_y, color="#0f0")

            actions.append(Deletion(last.end_x, last.end_y, deletion_length))
            plot.line(last.end_x, cur.start_y, cur.start_x, cur.start_y, color="#f00")

        elif cur.start_x < last.end_x and cur.start_y >= last.end_y:  # top left
            tmp_dot_y = YCoordOnLine(*last.coords, cur.start_x)
            insertion_length = cur.start_y - last.end_y
            duplication_length = last.end_x - cur.start_x
            duplication_height = last.end_y - tmp_dot_y

            actions.append(Insertion(cur.start_x, last.end_y, insertion_length))
            plot.line(cur.start_x, last.end_y, cur.start_x, cur.start_y, color="#0f0")

            actions.append(Duplication(cur.start_x, tmp_dot_y, duplication_length, duplication_height, line_index - 1))
            plot.poligon([
                (cur.start_x, tmp_dot_y),
                (cur.start_x, last.end_y),
                (last.end_x, last.end_y)
            ], color="#f0f")

        elif cur.start_x >= last.end_x and cur.start_y < last.end_y:  # bottom right
            deletion_length = cur.start_x - last.end_x
            translocation_length = last.end_y - cur.start_y

            actions.append(Deletion(last.end_x, last.end_y, deletion_length))
            plot.line(last.end_x, last.end_y, cur.start_x, last.end_y, color="#f00")

            actions.append(Translocation(last.end_x, last.end_y, translocation_length))
            plot.line(cur.start_x, last.end_y, cur.start_x, cur.start_y, color="#0ff")

        else:
            print("\nUnknown!!!\n")

        if cur.end_x >= last.end_x:
            last = cur

    large_actions = [action for action in actions if action.size >= settings["min_event_size"]]
    large_actions.sort(key=lambda action: -action.size)

    print("\nActions:", *actions, sep='\n')
    print("\nLarge_actions:", *large_actions, sep='\n')
    print()

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

    large_actions = [Pass()] + rotation_actions + large_actions

    with open(mkpath(output_folder, "history.txt"), 'w', encoding="utf-8") as history_file:
        for action in large_actions:

            if isinstance(action, Rotation):
                print("Rotation from {} (Query) to {} (Query)\n".format(
                    prtNum(int(lines[action.start_line].start_x)),
                    prtNum(int(lines[action.end_line].end_x))
                ), file=history_file)

            elif isinstance(action, Deletion):
                print("Deletion of {}-{} (Query) from {} (Ref)\n".format(
                    prtNum(int(action.start_x)),
                    prtNum(int(action.start_x + action.length)),
                    prtNum(int(action.start_y))
                ), file=history_file)

            elif isinstance(action, Insertion):
                print("Insertion of {}-{} (Ref) to {} (Query)\n".format(
                    prtNum(int(action.start_y)),
                    prtNum(int(action.start_y + action.height)),
                    prtNum(int(action.start_x))
                ), file=history_file)

            elif isinstance(action, Translocation):
                print("Translocation of {}-END (Query) from {} (Ref) to {} (Ref)\n".format(
                    prtNum(int(action.start_x)),
                    prtNum(int(action.start_y - action.height)),
                    prtNum(int(action.start_y))
                ), file=history_file)

            elif isinstance(action, Duplication):
                print("Duplication of {}-{} (Query) {}-{} (Ref)\n".format(
                    prtNum(int(action.start_x)),
                    prtNum(int(action.start_x + action.length)),
                    prtNum(int(action.start_y)),
                    prtNum(int(action.start_y + action.height))
                ), file=history_file)

    print(" {} images\n".format(len(large_actions)))

    # print("Large actions:", *large_actions, sep='\n')

    for action_index, action in enumerate(large_actions):

        if isinstance(action, Rotation):

            for i in range(action.start_line, action.end_line + 1):
                for j in range(len(lines[i].dots)):
                    lines[i].dots[j][1] -= (lines[i].dots[j][1] - action.rotation_center) * 2

            # for i in range(len(dots)):
            #     if start_query_pos <= dots[i][0] <= start_query_pos + length:
            #         dots[i][1] = dots[i][1] - (dots[i][1] - action.rotation_center) * 2

        elif isinstance(action, Insertion):
            for line in rotated_lines:
                new_dots = []
                for dot_x, dot_y in line.dots:
                    if dot_x > action.start_x:
                        dot_y -= action.height
                    if dot_x != action.start_x:
                        new_dots.append([dot_x, dot_y])

                line.dots = new_dots

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i].start_x >= action.start_x:
                    large_actions[i].start_y -= action.height

        elif isinstance(action, Deletion):
            for line in rotated_lines:
                for i in range(len(line.dots)):
                    if line.dots[i][0] >= action.start_x + action.length:
                        line.dots[i][0] -= action.length

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i].start_x >= action.start_x + action.length:
                    large_actions[i].start_x -= action.length

        elif isinstance(action, Duplication):
            new_dots = []
            for dot in rotated_lines[action.line_index].dots:
                if not (action.start_x <= dot[0] <= action.start_x + action.length):
                    new_dots.append(dot)
            rotated_lines[action.line_index].dots = new_dots

            for line in rotated_lines:
                for i in range(len(line.dots)):
                    if line.dots[i][0] >= action.start_x:
                        line.dots[i][1] -= action.height

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i].start_x >= action.start_x:
                    large_actions[i].start_y -= action.height

        elif isinstance(action, Translocation):
            for line in rotated_lines:
                new_dots = []
                for dot_x, dot_y in line.dots:
                    if dot_x > action.start_x:
                        dot_y += action.height
                    if dot_x != action.start_x:
                        new_dots.append([dot_x, dot_y])

                line.dots = new_dots

            for i in range(action_index + 1, len(large_actions)):
                if large_actions[i].start_x >= action.start_x:
                    large_actions[i].start_y += action.height

        elif isinstance(action, Pass):
            pass

        else:
            raise ValueError("History: Unknown action type")

        if isinstance(action, (Pass, Rotation)):
            # plot.scatter(dots, dotsize=settings["dotsize"], color="#00f")
            for line in lines:
                plot.scatter(line.dots, dotsize=settings["dotsize"], color="#00f")
        else:
            # Adjusting axes (bottom):
            bottom = float("inf")
            for line in rotated_lines:
                for dot_x, dot_y in line.dots:
                    bottom = min(bottom, dot_y)

            for line in rotated_lines:
                for i in range(len(line.dots)):
                    line.dots[i][1] -= bottom

            for i in range(action_index + 1, len(large_actions)):
                if hasattr(large_actions[i], "start_x") and hasattr(large_actions[i], "start_y") and \
                        large_actions[i].start_x >= action.start_x:
                    large_actions[i].start_y += bottom

            for line in rotated_lines:
                plot.scatter(line.dots, dotsize=settings["dotsize"], color="#00f")

        print("Saving large action #{}{}...\n".format(action_index, "" if isinstance(action, Pass) else " ({})".format(action.type)))
        plot.tight()
        plot.save(mkpath(
            output_folder,
            "history",
            "{}{}.png".format(
                str(action_index).zfill(3),
                "" if isinstance(action, Pass) else " ({})".format(action.type)
            )
        ))
        plot.clear()

    del plot


if __name__ == "__main__":
    analyze(query_genome_path, ref_genome_path, sam_file_path, show_plot, output_folder, SETTINGS)
