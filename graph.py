import matplotlib.pyplot as plt

FONT_SIZE = 6
DOT_SIZE = 10
DIRECTIONS8 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
DIRECTIONS4D = [(-1, -1), (-1, 1), (1, -1), (1, 1)]


def main(genome1_path, genome2_path):
    # with open(pairwise_path, 'r', encoding="utf-8") as pairwise_file:
    #     info = pairwise_file.readline()
    #     genome1 = pairwise_file.readline()
    #     alignment = pairwise_file.readline()
    #     genome2 = pairwise_file.readline()

    with open(genome1_path, 'r', encoding="utf-8") as file:
        genome1 = "".join(list(map(lambda line: line.strip(), file.readlines()[1:])))

    with open(genome2_path, 'r', encoding="utf-8") as file:
        genome2 = "".join(list(map(lambda line: line.strip(), file.readlines()[1:])))

    length = max(len(genome1), len(genome2))
    genome1 = genome1 + "-" * (length - len(genome1))
    genome2 = genome2 + "-" * (length - len(genome2))

    matrix = [[False] * length for _ in range(length)]

    # dots = []

    for i in range(length):
        for j in range(length):
            if genome1[i] == genome2[j]:
                matrix[i][j] = True

    for i in range(1, length - 1):
        for j in range(1, length - 1):
            if matrix[i][j]:
                for d_x, d_y in DIRECTIONS8:
                    if matrix[i + d_y][j + d_x]:
                        break
                else:
                    matrix[i][j] = False

    for i in range(1, length - 1):
        for j in range(1, length - 1):
            if matrix[i][j]:
                for d_x, d_y in DIRECTIONS4D:
                    if matrix[i + d_y][j + d_x]:
                        break
                else:
                    matrix[i][j] = False

    X = 10

    def go(i, j, count_first, count_second, length, direction):
        pass

    for iteration in range(3):
        for i in range(X + 1, length - 1 - X):
            for j in range(X + 1, length - 1 - X):
                if matrix[i][j]:

                    first = 0
                    second = 0
                    for x in range(X):
                        first += int(matrix[i + x][j + x] or matrix[i - x][j - x])
                        second += int(matrix[i - x][j + x] or matrix[i + x][j - x])

                    if first != X and second != X:
                        matrix[i][j] = False

    plt.rc('font', size=FONT_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=FONT_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=8)             # fontsize of the x and y labels
    plt.rc('xtick', labelsize=FONT_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=FONT_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=FONT_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=FONT_SIZE)   # fontsize of the figure title

    plt.xticks(rotation=-90)

    plt.imshow(matrix[::-1])

    # plt.scatter(*zip(*dots), s=DOT_SIZE).set_label('')
    plt.legend(loc="center left")
    plt.show()


# ---SETTINGS--- #

# pairwise_path = "BWA/small/deletion/bwa_output_pairwise.txt"
# pairwise_path = "BWA/large/bwa_output_pairwise.txt"

genome1_path = "samples/small/source.fasta"
genome2_path = "samples/small/inversion.fasta"

# ---SETTINGS--- #

if __name__ == "__main__":
    main(genome1_path, genome2_path)
