from sys import argv

from time import time

import bz2


def main():
    """
    """
    start = time()

    fname = argv[1]

    if '-bz2' in argv:
        count = count_line_bz2(fname)
    else:
        count = count_line(fname)

    print("count done in {0} s".format(time() - start))
    print("number of lines {0}".format(count))


def count_line(fname):
    """
    """
    f_bz2 = open(fname)

    count = 0

    for line in f_bz2:
        count += 1

    return count


def count_line_bz2(fname):
    """
    """
    f_bz2 = bz2.BZ2File(fname)

    count = 0

    for line in f_bz2:
        count += 1

    return count


if __name__ == '__main__':
    main()
