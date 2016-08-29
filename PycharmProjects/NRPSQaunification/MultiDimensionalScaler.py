from sklearn import manifold
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
import numpy as np

def create_mds(data_matrix, name_list):
    mds = manifold.MDS(n_components=2, max_iter=100, n_init=1)
    pos = mds.fit_transform(data_matrix)
    fig = plt.figure(1)
    ax = plt.axes([0., 0., 1., 1.])
    plt.scatter(pos[:, 0], pos[:, 1], s=20, c='g')
    for name, x, y in zip(name_list, pos[:, 0], pos[:, 1]):
        plt.annotate(name, xy = (x, y), xytext=(-20, 20), textcoords='offset points', ha='right', va='bottom')
    plt.legend(('MDS'), loc='best')
    plt.show()
