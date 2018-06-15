import seaborn as sns
import pandas as pd

from glob import glob

import pylab as plt

from collections import Counter

import numpy as np


################ VARIABLE ################################
PROJECT_NAME = '10p7reads'
NO_FILTER = "NOCORRECTIONs"
PATH_DATA = '/home/oliver/data/ATACtest/'
##########################################################


def main():
    """
    """
    fig, ax = plt.subplots(figsize=(12, 12))

    fname = glob('{0}/*{1}.success.log'.format(PATH_DATA, PROJECT_NAME))[0]
    data_success = load_files(fname)

    value_success = list(data_success.values())

    fname = glob('{0}/*{1}.fail.log'.format(PATH_DATA, PROJECT_NAME))[0]
    data_fail = load_files(fname)

    value_fail = list(data_fail.values())


    fname = glob('{0}/*{1}.success.log'.format(PATH_DATA, NO_FILTER))[0]
    data_ref = load_files(fname)

    value_ref = list(data_ref.values())

    frame_success = pd.DataFrame({'type':'success', "value": value_success})
    frame_fail = pd.DataFrame({'type':'fail', "value": value_fail})
    frame_ref = pd.DataFrame({'type':'ref', "value": value_ref})
    # frame = pd.concat([frame_success, frame_fail])

    sns.kdeplot(np.log10(1.0 + frame_success['value']), shade=True, ax=ax, label='success tag')
    sns.kdeplot(np.log10(1.0 + frame_fail['value']), shade=True, color='r', ax=ax, label='failed tag')
    # sns.kdeplot(np.log10(1.0 + frame_ref['value']), shade=True, color='y', ax=ax, label='No filter')

    ax.set_xlabel('log10(Number of tags)')
    ax.set_ylabel('Density')
    ax.set_title('Density plots of the number of reads per group of tags')

    for read in data_fail.most_common(n=20):
        print("tag: {0} ".format(read))

    plt.show()


def load_files(fname, sep='\t'):
    """
    """
    dic = Counter()

    for line in open(fname):
        line = line.strip('\n').split()
        dic[line[0].strip(':')] = int(line[1])

    return dic


if __name__ == '__main__':
    main()
