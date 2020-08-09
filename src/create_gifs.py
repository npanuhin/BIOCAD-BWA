from subprocess import Popen, PIPE
import os


def createGif(path_in, path_out, delay=20, resize=0.3):
    s = Popen(
        "convert -delay {} -loop 0 {} -resize {}% {}".format(delay, path_in, resize * 100, path_out),
        shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()

    print(s[0])
    print(s[1])


def main():
    if not os.path.isdir("history"):
        os.mkdir("history")

    os.chdir("history")

    print("Now in {}".format(os.getcwd()))

    createGif("*.png", "../history.gif")

    print("-" * 30 + "\n")

    os.chdir("../")


os.chdir("../tests/")

for foldername in os.listdir():
    if not os.path.isdir(foldername):
        continue

    os.chdir(foldername)

    if foldername.strip("/").strip("\\") == "small":
        for foldername2 in os.listdir():
            os.chdir(foldername2)
            main()
            os.chdir("../")
    else:
        main()

    os.chdir("../")
