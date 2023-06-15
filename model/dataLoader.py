import pandas as pd
import numpy as np

# min-max normalization
def MinMax(x):
    '''

    :param x:
    :return: min-max normalization x

    '''
    x = (x - np.min(x)) / (np.max(x) - np.min(x))
    return x


def dataLoad(filepath):
    '''

    :param filepath path of expression matrix
    :return: preprocess expression matrix and gene keys
    '''

    dat = pd.read_csv(filepath)
    gene = dat.iloc[:, 0].to_list()
    exp = dat.iloc[:, 1:dat.shape[1]].transpose()

    scale_exp = exp.apply(lambda x: MinMax(x), axis=1)

    return scale_exp, gene

