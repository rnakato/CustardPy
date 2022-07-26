import os
import sys
import numpy as np
import scipy.stats as sp
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from custardpy.InsulationScore import *
from custardpy.DirectionalityIndex import *
from custardpy.loadData import loadDenseMatrix

def getNonZeroMatrix(A, lim_pzero):
    A = A.fillna(0)
    pzero_row = A[A>0].count(axis=0)/A.shape[0]
    index = pzero_row[pzero_row > lim_pzero].index
    pzero_col = A[A>0].count(axis=1)/A.shape[1]
    columns = pzero_col[pzero_col > lim_pzero].index

    A = A[index]
    A = A.loc[columns]

    return A

class JuicerMatrix:
    def __init__(self, norm, rawmatrix, res, *, eigenfile=""):
        self.res = res
        if os.path.exists(rawmatrix):
            self.raw = loadDenseMatrix(rawmatrix)
        else:
            print ("Error: " + rawmatrix + " does not exist.")
            sys.exit()
        if os.path.exists(eigenfile):
            self.eigen = np.loadtxt(eigenfile)
            eigen_sortindex = np.argsort(self.eigen)
            eigen_sort = self.eigen[eigen_sortindex]
            self.eigen_sortindex = eigen_sortindex[~np.isnan(eigen_sort)]
#            eigen_sort = eigen[eigen_sortindex]
        else:
            self.eigen = np.zeros(self.raw.shape[0])
        if norm == "RPM":
            self.raw = self.raw * 10000000 / np.nansum(self.raw)
        self.InsulationScore = MultiInsulationScore(self.getmatrix().values,
                                                    1000000, 100000, self.res)

    def getmatrix(self, *, isNonZero=False, sortEigen=False):
        if isNonZero == True:
            return getNonZeroMatrix(self.raw, 0)
        elif sortEigen == True:
            m = self.raw.iloc[self.eigen_sortindex,:]
            m = m.iloc[:,self.eigen_sortindex]
            return m
        else:
            return self.raw

    def getlog(self, *, isNonZero=False, sortEigen=False):
        mat = self.getmatrix(isNonZero=isNonZero, sortEigen=sortEigen)
        logmat = mat.apply(np.log1p)
        return logmat

    def getZscore(self):
        logmat = self.getlog()
        zmat = pd.DataFrame(sp.stats.zscore(logmat, axis=1), index=logmat.index, columns=logmat.index)
        return zmat

#    def getPearson(self, *, isOE=False):
#        oe = self.getlog(isOE=isOE)
#        ccmat = np.corrcoef(oe)
#        ccmat[np.isnan(ccmat)] = 0
#        ccmat = pd.DataFrame(ccmat, index=oe.index, columns=oe.index)
#        return ccmat

    def getEigen(self, *, sortEigen=False):
        if sortEigen == True:
            return self.eigen[self.eigen_sortindex]
        else:
            return self.eigen
#        from sklearn.decomposition import PCA
#        pca = PCA()
#        ccmat = self.getPearson()
#        index = np.isnan(ccmat).all(axis=1)
#        transformed = pca.fit_transform(ccmat)
#        pc1 = transformed[:, 0]
#        pc1[index] = np.nan
#        return transformed[:, 0]

    def getInsulationScore(self, *, distance=500000):
        return self.InsulationScore.getInsulationScore(distance=distance)

    def getMultiInsulationScore(self):
        return self.InsulationScore.getMultiInsulationScore()

    def getDirectionalityIndex(self, *, distance=1000000):
        return calcDI(self.raw.values, self.res, distance=distance)

def ExtractMatrix(mat,s,e):
    if e==-1:
        return mat.values[s:,s:]
    else:
        return mat.values[s:e,s:e]



def ExtractMatrixIndex(mat,index1,index2):
    mat = mat[index1,:]
    mat = mat[:,index2]
    return mat

def ExtractTopOfPC1(mat, pc1, nbin):
    sortedindex = np.argsort(pc1)
    sortmat = mat[sortedindex,:]
    sortmat = sortmat[:,sortedindex]
    sortmat = np.concatenate((sortmat[:nbin,:] , sortmat[-nbin:,:]), axis=0)
    sortmat = np.concatenate((sortmat[:,:nbin] , sortmat[:,-nbin:]), axis=1)
    return sortmat

def getCompartment(mat, pc1):
    indexA = pc1 > 1
    indexB = pc1 < -1
    A = mat[indexA]
    A = A[:,indexA]
    B = mat[indexB]
    B = B[:,indexB]
    return A, B

def getCompartment_inter(mat, pc1_odd, pc1_even, nbin):
    sortedindex_odd = np.argsort(pc1_odd)
    sortedindex_even = np.argsort(pc1_even)
    sortmat = mat[sortedindex_odd,:]
    sortmat = sortmat[:,sortedindex_even]
    A = sortmat[:nbin,:nbin]
    B = sortmat[-nbin:,-nbin:]
#    indexA_odd = pc1_odd > 1
#    indexB_odd = pc1_odd < -1
 #   indexA_even = pc1_even > 1
  #  indexB_even = pc1_even < -1
  #  A = mat[indexA_odd]
  #  A = A[:,indexA_even]
  #  B = mat[indexB_odd]
   # B = B[:,indexB_even]
    return A, B
