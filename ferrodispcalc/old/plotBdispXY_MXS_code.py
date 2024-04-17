#!/usr/bin/python
import numpy as np
from math import pi, cos, sqrt
import matplotlib
import sys
matplotlib.use('Agg')
import matplotlib.pylab as plt
from matplotlib import cm
import os
import subprocess

#Usuage: python plotDumpLocalP.py dump.localP 10
#Plot the dipole configurations every 10th step

# Supercell size: NX*NY*NZ
NX = int(sys.argv[2])
NY = int(sys.argv[3])
NZ = int(sys.argv[4])
n_A = NX * NY * NZ


def calAngle(dx, dy):
    pp = np.sqrt(dx * dx + dy * dy)
    angle = np.arccos(dx / pp) / pi * 180.0
    index = np.where(dy < 0.0)
    angle[index] = 360.0 - angle[index]
    return angle


def plotXZ(ii, dxx, dyy):
    dx0 = dxx[ii, :, :]
    dy0 = dyy[ii, :, :]
    fig, ax = plt.subplots()
    angle = calAngle(dx0.flatten(), dy0.flatten())
    profile = angle.reshape(dx0.shape)
    ax.quiver(dx0, dy0)
    sc = ax.imshow(profile,
                   cmap=cm.hsv,
                   vmax=360,
                   vmin=0,
                   aspect=1.0,
                   origin='lower')
    plt.savefig("%s_XY/Profile_Layer%2.2d_xy.png" %
                (sys.argv[1], ii),
                bbox_inches='tight',
                dpi=300)
    plt.close()

os.system("rm -rf %s_XY" % sys.argv[1])
os.mkdir("%s_XY" % sys.argv[1])

skipN = 5 + n_A
data = np.loadtxt("./avg_disp", skiprows=skipN, usecols=(0,1,2))

dx = data[:, 0]
dy = data[:, 1]
dz = data[:, 2]

dx = dx.reshape((NZ, NY, NX))
dy = dy.reshape((NZ, NY, NX))
dz = dz.reshape((NZ, NY, NX))

for jj in range(0, NZ):
    print("Layer:", jj)
    plotXZ(jj, dx, dy)
