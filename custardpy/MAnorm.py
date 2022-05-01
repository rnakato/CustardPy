import numpy as np
import scipy.stats as sp
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import linear_model
#from module.HiCmodule import *

class MAdata: # log, asinhを仮定
    def __init__(self, mat1, mat2):
        self.mat1 = mat1
        self.mat2 = mat2
        mat1 = np.triu(mat1)
        mat2 = np.triu(mat2)
        self.M = mat2 - mat1
        self.A = (mat1 + mat2)/2

        index = np.where((mat1 > 0.1) & (mat2 > 0.1))
        mat1 = mat1[index]
        mat2 = mat2[index]
        index = np.where((mat1-mat2 > -2.5) & (mat1-mat2 < 2.5))
        mat1 = mat1[index]
        mat2 = mat2[index]
        self.Mlim = mat2 - mat1
        self.Alim = (mat1 + mat2)/2
        self.clf = linear_model.LinearRegression(fit_intercept=True, normalize=True, copy_X=True, n_jobs=-1)
        self.clf.fit(self.Alim.reshape(-1, 1), self.Mlim)

        self.Aadjust = self.Alim
        self.Madjust = self.Mlim - self.clf.coef_ * self.Alim - self.clf.intercept_

    def plot(self, ax, label, *, isMAnorm=False, isLineFitted=False):
        from matplotlib import lines
#        if plottype == "trim":
 #           ax.scatter(self.Alim, self.Mlim, alpha = 0.3)
        if isMAnorm == True:
            ax.scatter(self.Aadjust, self.Madjust, alpha = 0.3)
        else:
            ax.scatter(self.A, self.M, alpha = 0.3)
        if isLineFitted == True:
            line = lines.Line2D([0, 8], [self.clf.intercept_, self.clf.intercept_ + 8*self.clf.coef_], color='r', linewidth=2)
            ax.add_line(line)
        else:
            ax.hlines(y=0, xmin=0, xmax=8, colors='r', linewidths=2)
        ax.set_title(label)

    def plotMAnorm(self):
        fig = plt.figure(figsize=(14, 3))
        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2)
        ax3 = fig.add_subplot(1,3,3)
        self.plot(ax1, "raw")
        self.plot(ax2, "fit", isLineFitted=True)
        self.plot(ax3, "MAnorm", isMAnorm=True)

    def getM(self):
        return self.mat2 - self.mat1

    def getMadjust(self):
        M = self.getM()
        A = (self.mat1 + self.mat2)/2
        Madjust = M - self.clf.coef_*A - self.clf.intercept_
        return Madjust

    def plotRatioMatrix(self):
        cm = generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f'])
        fig = plt.figure(figsize=(7, 7))
        ax1 = fig.add_subplot(2,2,1)
        ax2 = fig.add_subplot(2,2,2)
        ax3 = fig.add_subplot(2,2,3)
        ax4 = fig.add_subplot(2,2,4)
        self.plot(ax1, "fit", isLineFitted=True)
        self.plot(ax2, "MAnorm", isMAnorm=True)
        ax3.imshow(self.getM(), clim=(-3, 3), cmap=cm)
        ax4.imshow(self.getMadjust(), clim=(-3, 3), cmap=cm)
