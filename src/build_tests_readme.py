import os


def buildTests():
    os.chdir("tests")

    with open("README.md", 'w', encoding="utf-8") as readme:
        print("<h1 align=\"center\">Tests</h1>\n", file=readme)

        test_folders = [foldername for foldername in os.listdir() if os.path.isdir(foldername) and foldername != "small"]

        headers = ["|", "BLAST", "My plot", "History GIF", "History text output"]
        table = []

        for test_folder in test_folders:
            # os.chdir(test_folder)

            blast_link = "https://raw.githubusercontent.com/npanuhin/BIOCAD_BWA/master/tests/{}/BLAST.png".format(test_folder)
            plot_link = "https://raw.githubusercontent.com/npanuhin/BIOCAD_BWA/master/tests/{}/sam_analyze.png".format(test_folder)
            histor_gif_link = "https://raw.githubusercontent.com/npanuhin/BIOCAD_BWA/master/tests/{}/history.gif".format(test_folder)
            history_text_link = "https://github.com/npanuhin/BIOCAD_BWA/blob/master/tests/{}/history.txt".format(test_folder)

            table.append([
                test_folder,
                "[BLAST]({} \"View image\")".format(blast_link.replace(' ', "%20")),
                "[My plot]({} \"View image\")".format(plot_link.replace(' ', "%20")),
                "[History GIF]({} \"View GIF\")".format(histor_gif_link.replace(' ', "%20")),
                "[History text output]({} \"View JSON file\")".format(history_text_link.replace(' ', "%20"))
            ])

            # os.chdir("../")

        print("|".join(headers), file=readme)
        print("|".join([":-:"] * len(headers)), file=readme)
        for line in table:
            print("|".join(line), file=readme)

    os.chdir("../")


def main(root_path):
    start_path = os.getcwd()
    os.chdir(root_path)
    buildTests()
    os.chdir(start_path)


if __name__ == "__main__":
    main("../")
