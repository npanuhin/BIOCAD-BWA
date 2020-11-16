from os import walk as os_walk
from os.path import normpath as os_normpath, join as os_join
from shutil import rmtree


def mkpath(*paths):
    return os_normpath(os_join(*paths))


def removePythonCache(root="../"):
    for cur_path, folders, files in os_walk(root):
        for folder in folders:
            if folder == "__pycache__":
                rmtree(mkpath(cur_path, folder))


if __name__ == "__main__":
    removePythonCache()
