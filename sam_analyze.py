from string import ascii_uppercase


def prettifyNumber(num):
    return "{:,}".format(num)

def main(sam_file_path):

    bitwise_flags = [
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

    with open(sam_file_path, 'r', encoding="utf-8") as sam_file:

        line_index = -1
        sam_info = []

        sam_data = []

        for line in sam_file:
            line_index += 1

            if line.startswith('@'):
                sam_info.append(line)
                continue

            line = line.split()
            mapQ = int(line[4])
            quality = round((10 ** (mapQ / -10)) * 100, 5)
            position = int(line[3])
            flags_bit = int(line[1])
            flags = []

            for i in range(len(bitwise_flags) - 1, -1, -1):
                cur_flag, flags_bit = divmod(flags_bit, 2 ** i)
                if cur_flag:
                    flags.append(i)

            flags.sort()

            rid = line[9]
            rid_size = len(rid)

            actions = line[5]

            if rid_size <= 10000:
                continue

            sam_data.append([position, flags, quality, rid_size, actions])

        sam_data.sort(key=lambda item: item[2])

        for position, flags, quality, rid_size, actions in sam_data:

            print("{} -> {} [Q:{}]".format(prettifyNumber(position), prettifyNumber(position + rid_size), quality))

            for flag in flags:
                if flag == 11:  # Currently flag №11 is disabled
                    continue
                print(bitwise_flags[flag])

            for i in range(len(actions)):
                if actions[i] in ascii_uppercase:
                    print("-{}|".format(actions[i]), end="")
                else:
                    print(actions[i], end="")
            print()

            print()


# ---SETTINGS--- #

sam_file_path = "BWA/large3/bwa_output.sam"

# ---SETTINGS--- #

if __name__ == "__main__":
    main(sam_file_path)
