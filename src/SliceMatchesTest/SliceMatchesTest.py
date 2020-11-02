from Plot import Plot


# ---SETTINS--- #

slice_matches_file_path = "test_input.txt"

# ---SETTINS--- #

def main(slice_matches_file_path):

    with open(slice_matches_file_path, 'r', encoding="utf-8") as input_file:
        _, query_name = map(str.strip, input_file.readline().split(':', 1))
        _, target_name = map(str.strip, input_file.readline().split(':', 1))

        query_name, target_name = query_name.strip(), target_name.strip()

        print("Query: {}".format(query_name))
        print("Reference: {}".format(target_name))

        plot = Plot("SliceMatchesTest", 8, grid_size=None, figsize=(10, 7), nameX=query_name, nameY=target_name)

        for line in input_file:
            block, positions = map(str.strip, line.split(':'))

            query_positions, target_positions = map(str.strip, positions.split(';'))

            query_positions = list(map(int, query_positions.split(',')))
            target_positions = list(map(int, target_positions.split(',')))

            block_size = len(block)

            # print(block, query_positions, target_positions)

            for query_pos in query_positions:
                for target_pos in target_positions:
                    plot.line(query_pos, target_pos, query_pos + block_size, target_pos + block_size, color="#000")

    plot.show()
    del plot


if __name__ == "__main__":
    main(slice_matches_file_path)
