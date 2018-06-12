from sys import argv

from time import time

import bz2


def main():
    """
    """
    start = time()

    f_bz2 = open(argv[1])

    count = 0

    for line in f_bz2:
        count += 1

    print("count done in {0} s".format(time() - start))
    print("number of lines {0}".format(count))


if __name__ == '__main__':
    main()
