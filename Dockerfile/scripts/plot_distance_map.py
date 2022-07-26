#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import numba
# --------------------------------------------------------------------------------------------------
# Set command line arguments
argv = sys.argv
argc = len(argv)
if (argc != 5):
    print("Usage: python " + argv[0] + " FILE_READ_XYZ N FRAME ODIR")
    exit()
FILE_READ = argv[1]
N = int(argv[2])
FRAME = int(argv[3])
# --------------------------------------------------------------------------------------------------
DIR = argv[4]
os.makedirs(DIR, exist_ok=True)
# --------------------------------------------------------------------------------------------------


@numba.jit
def Calc_Distance(xi, yi, zi, xj, yj, zj):
    distance = np.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)
    return distance
# --------------------------------------------------------------------------------------------------


@numba.jit
def Calc_Distance_Map(Rx, Ry, Rz):
    DM = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            DM[i, j] = Calc_Distance(Rx[i], Ry[i], Rz[i], Rx[j], Ry[j], Rz[j])
            DM[j, i] = DM[i, j]
    return DM
# --------------------------------------------------------------------------------------------------


def main():
    Rx = np.zeros(N)
    Ry = np.zeros(N)
    Rz = np.zeros(N)
    # ----------------------------------------------------------------------------------------------
    fp = open(FILE_READ, "r")
    # ----------------------------------------------------------------------------------------------
    for frame in range(FRAME + 1):
        line = fp.readline()
        line = fp.readline()
        # ------------------------------------------------------------------------------------------
        for n in range(N):
            line = fp.readline()
            items = line.split("\t")
            if(len(items)<3):
                continue
            Rx[n] = float(items[1])
            Ry[n] = float(items[2])
            Rz[n] = float(items[3])
        # ------------------------------------------------------------------------------------------
        # if frame == 662 or frame == 665 or frame == 668 or frame == 671:
        DM = Calc_Distance_Map(Rx, Ry, Rz)
        # ------------------------------------------------------------------------------------------
        plt.rcParams["font.family"] = "Arial"
        plt.rcParams["font.size"] = 30
        # ------------------------------------------------------------------------------------------
#        FILE_FIG = DIR + "/frame_{0:04d}.png".format(frame)
#        FILE_FIG = DIR + "/frame_{0:04d}.svg".format(frame)
        FILE_FIG = DIR + "/frame_{0:04d}.jpg".format(frame)
        plt.figure(figsize=(6, 5))
        plt.imshow(DM, cmap="coolwarm_r", clim=(0, 4))
        plt.colorbar(ticks=[0, 1, 2, 3, 4], shrink=0.6)
        # plt.imshow(DM, cmap="coolwarm_r", clim=(0, 8))
        # plt.colorbar(ticks=[0, 2, 4, 6, 8], shrink=0.6)
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(FILE_FIG)
        plt.close()
        # ------------------------------------------------------------------------------------------
    fp.close()
# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
