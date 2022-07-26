import numpy as np
import matplotlib.pyplot as plt

def getlinalg(gradient, wid):
    data = []
    for i in np.arange(gradient.size - wid):
        array = gradient[i:i+wid]
        data.append(np.linalg.norm(array))
    return data

def getMaxGradient(arr, val, order_argrelmax, wid, posi):
    from scipy.signal import argrelmax
    index = argrelmax(arr, order=order_argrelmax)
    index = index[0][np.where(arr[index[0]] > val)]
    index += wid # linalgのindexは元配列よりwidだけずれているため
    index = index[np.where((index < posi - 5) | (index > posi + 5))] # 対角線周辺(±5bin)の勾配は無視
    return index

def imshowWithvLine(A, posi, index, title, cm):
    fig = plt.figure(figsize=(6, 6))
    ax1 = fig.add_subplot(1,1,1)
    ax1.imshow(A, clim=(-2, 2), cmap=cm)
    ax1.set_title(title)
    plt.hlines(y=posi, xmin=0, xmax=190, colors='black', linewidths=1)
    for x in index:
        plt.vlines(x=x, ymin=0, ymax=190, colors='green', linewidths=1, linestyles='dashed')

def getIndexMatrix(linalg, segment_len, limit_val, order_argrelmax, wid):
    indexMatrix = []
    length = int(linalg.shape[0]/segment_len)
    for i in range(length):
        lin = linalg[i*segment_len:(i+1)*segment_len].mean(axis=0)
        index_merged = getMaxGradient(lin, limit_val, order_argrelmax, wid, (i+0.5)*segment_len)
        indexMatrix.append(index_merged)
        
    return indexMatrix

def imshowWithBoundary(A, indexMatrix, segment_len, cm):
    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(1,1,1)
    ax1.imshow(A, clim=(-2, 2), cmap=cm)
    for x, array in enumerate(indexMatrix):
        for y in array:
            plt.vlines(x=y, ymin=x*segment_len, ymax=(x+1)*segment_len, colors='green', linewidths=1)
    for x, array in enumerate(indexMatrix):
        for y in array:
            plt.hlines(y=y, xmin=x*segment_len, xmax=(x+1)*segment_len, colors='green', linewidths=1)

def getPeakfromIndexMatrix(indexMatrix, order):
    import collections
    flatten = []
    for array in indexMatrix:
        flatten.extend(array)

    dict = collections.Counter(flatten)
    lists = sorted(dict.items())
    x, y = zip(*lists)
    array = np.zeros(max(x)+1)
    for a, b in zip(x,y):
        array[a] = b

    from scipy.signal import argrelmax
    index = argrelmax(array, order=order)
    index = index[0]
    
    fig = plt.figure(figsize=(12, 4))
    plt.plot(array)
    for x in index:
        plt.vlines(x=x, ymin=0, ymax=array.max(), colors='green', linewidths=1, linestyles='dashed')
    plt.show()
    
    return index
