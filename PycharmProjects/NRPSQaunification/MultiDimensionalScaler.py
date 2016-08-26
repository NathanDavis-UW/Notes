from sklearn import manifold
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
import numpy as np

def create_mds(data_matrix, name_list):
    data_matrix = np.matrix(data_matrix)
    seed = np.random.RandomState(seed=3)
    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state = seed, dissimilarity="precomputed",
                  n_jobs=1)
    pos = mds.fit(data_matrix).embedding_
    #pos *= np.sqrt((data_matrix ** 2).sum()) / np.sqrt((pos ** 2).sum())
    clf = PCA(n_components=2)

    pos = clf.fit_transform(pos)
    fig = plt.figure(1)
    ax = plt.axes([0., 0., 1., 1.])
    plt.scatter(pos[: 0], pos[:, 1], s=20, c='g')
    plt.legend(('MDS'), loc='best')
    plt.show()
