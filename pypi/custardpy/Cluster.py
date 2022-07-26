import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def KMeansPlot(matrix, title, ncluster):
    import matplotlib.cm
    from sklearn.cluster import MiniBatchKMeans
    model = MiniBatchKMeans(random_state=0, n_clusters=ncluster, max_iter=10000)
    kmeans = model.fit_predict(matrix)
    color = matplotlib.cm.brg(np.linspace(0,1, np.max(kmeans) - np.min(kmeans)+1))

    for i in range(np.min(kmeans), np.max(kmeans)+1):
        plt.plot(matrix[kmeans == i][:,0],    matrix[kmeans == i][:,1], ".", color=color[i])
        plt.text(matrix[kmeans == i][:,0][0], matrix[kmeans == i][:,1][0], str(i+1), color="black", size=16)
    plt.title(title, size=16)

    return kmeans

# DBSCANでクラスタリング
def DBSCANPlot(matrix, title, eps):
    from sklearn.cluster import DBSCAN
    import matplotlib.cm
    model = DBSCAN(eps=eps)
    dbscan = model.fit_predict(matrix)
    color = matplotlib.cm.brg(np.linspace(0,1,np.max(dbscan) - np.min(dbscan)+1))

    for i in range(np.min(dbscan), np.max(dbscan)+1):
        plt.plot(matrix[dbscan == i][:,0],    matrix[dbscan == i][:,1], ".", color=color[i])
        plt.text(matrix[dbscan == i][:,0][0], matrix[dbscan == i][:,1][0], str(i+1), color="black", size=16)
    plt.title(title, size=16)
    return dbscan

def getSumMatrix(A, boundary):
    submatrix = np.vsplit(A, boundary)
    for i, mat in enumerate(submatrix):
        ms = np.hsplit(mat, boundary)
        for j, m in enumerate(ms):
            if j==0:
                data = m.mean()
            else:
                data = np.r_[data, m.mean()]
        if i==0:
            gs = data
        else:
            gs = np.c_[gs, data]
    return gs

def get_ellipse_coords(a=0.0, b=0.0, x=0.0, y=0.0, angle=0.0, k=2):
    """ Draws an ellipse using (360*k + 1) discrete points
    k = 1 means 361 points (degree by degree)
    a = major axis distance,
    b = minor axis distance,
    x = offset along the x-axis
    y = offset along the y-axis
    angle = clockwise rotation [in degrees] of the ellipse;
        * angle=0  : the ellipse is aligned with the positive x-axis
        * angle=30 : rotated 30 degrees clockwise from positive x-axis

    this function is obtained from : http://scipy-central.org/item/23/2/plot-an-ellipse
    """
    pts = np.zeros((int(180*k+1), 2))

    beta = -angle * np.pi/180.0
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)
    alpha = -np.radians(np.r_[0.:180.:1j*(180*k+1)])

    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)

    pts[:, 0] = x + (a * cos_alpha * cos_beta - b * sin_alpha * sin_beta)
    pts[:, 1] = y + (a * cos_alpha * sin_beta + b * sin_alpha * cos_beta)

    return pts

def restore_mat(mat, ref_data, columnname):
    a = pd.concat([ref_data, mat], axis=1, join='outer')
    a = a.iloc[:, ref_data.shape[1]:a.shape[1]]
    return a[columnname].unstack()

def plotArc(s, e):
    rad = (e - s)/2
    center = (s + e)/2
    pts = get_ellipse_coords(a=rad, b=1.0, x=center)
    plt.plot(pts[:,0], pts[:,1])

def plotVsegmentArc(vsegment, s, e, xstart, resolution):
    def getBed(index, xstart, resolution):
        s = index[0]/resolution - xstart
        e = index[-1]/resolution + 1 - xstart
        common = int((s+e)/2)
        return s,e,common

    if(s>e): return
    sindex = vsegment[s]
    eindex = vsegment[e]
    s1, e1, c1 = getBed(sindex, xstart, resolution)
    s2, e2, c2 = getBed(eindex, xstart, resolution)
    plotArc(c1, c2)
#    plt.axhline(y=0, xmin=s1, xmax=e1)
#    plt.axhline(y=0, xmin=s2, xmax=e2)

def plotVArc(s, e, xstart, resolution):
    def getBed(index, xstart, resolution):
        s = index/resolution - xstart
        e = index/resolution + 1 - xstart
        common = int((s+e)/2)
        return s,e,common

    if(s>e): return
    s1, e1, c1 = getBed(s, xstart, resolution)
    s2, e2, c2 = getBed(e, xstart, resolution)
    plotArc(c1, c2)
    plt.axhline(y=0, xmin=s1, xmax=e1)
    plt.axhline(y=0, xmin=s2, xmax=e2)

def getHead(labels):
    nlabels = len(labels)
    a = []
    for i in range(nlabels):
        for j in range(i+1,nlabels):
            a.append(labels[i] + "-" + labels[j])
    return a

def get_corr_allcluster(sub3d, ncluster, kmeans, labels):
    def getmat(df, nlabels):
        corr_mat = df.corr(method='spearman')
        corr_mat = corr_mat.values
        a = []
        for i in range(nlabels-1):
            for x in corr_mat[i,i+1:nlabels]: a.append(x)
        return a

    nlabels = len(labels)
    for cl in range(ncluster):
        df = pd.DataFrame(sub3d[kmeans==cl])
        a = getmat(df, nlabels)
        if cl==0:
            mat = a
        else:
            mat = np.c_[mat, a]

    df = pd.DataFrame(sub3d)
    a = getmat(df, nlabels)
    mat = np.c_[mat, a]

    corr_allcluster = pd.DataFrame(mat.T, columns=getHead(labels))
    corr_allcluster = corr_allcluster.rename(index={ncluster: 'All'})
    return corr_allcluster


def draw_samples_whole(matrix, labels, cm):
    nsample = matrix.shape[0]
    plt.figure(figsize=(16, 6))
    for i in range(nsample):
        ax = plt.subplot(1, nsample+1, i+1)
        plt.imshow(matrix[i], clim=(-2, 2), cmap=cm)
        ax.set_title(labels[i])
    plt.tight_layout()
    plt.show()

def addzero_to_3dmatrix(mat, difflength, resolution):
    diff = difflength / resolution
    lim_pzero = 0.1
    mat[np.isnan(mat)] = 0
    index_zero = np.sum(np.sum(mat, axis=0)>0, axis=1)/mat.shape[1] < lim_pzero
    mat[:,index_zero] = 0
    mat[:,:,index_zero] = 0
    for i in range(mat.shape[1]):
        for j in range(mat.shape[1]):
            if j < i: mat[:,i,j] = 0
            if j > i + diff: mat[:,i,j] = 0

def make_refdata(ref_matrix, ref, labels, resolution, difflim, cm):
    import copy
    matrix = copy.deepcopy(ref_matrix)
    # 0が多い行・列を0に
    # 10M以上離れた領域も0に
    addzero_to_3dmatrix(matrix, 10000000, resolution)

    draw_samples_whole(matrix, labels, cm)

    ref_data = matrix.reshape(ref_matrix.shape[0], ref_matrix.shape[1]*ref_matrix.shape[2]).T
    ref_data = pd.DataFrame(ref_data, index=ref.index, columns=labels)
    ref_data.index.names = ['position1', 'position2']
    return ref_data

def make_sub3d(_ref_data):
    # 0を１つでも含む行を削除
    mat = _ref_data
    nonzero = np.logical_not((mat == 0).any(axis=1))
    mat = mat[nonzero]
#    plt.imshow(restore_mat(sub3d, _ref_data, "Rad21").iloc[550:650,550:650], clim=(-2, 2), cmap=cm)
#    plt.show()
    return mat
