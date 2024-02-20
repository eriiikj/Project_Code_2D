#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:29:59 2024

@author: Erik Jacobsson
"""

from matplotlib import pyplot as plt
import numpy as np

x = np.array([1321.4, 598.6, 580.6, 563.8, 548.6, 535.4, 524.5, 516.2, 511,
509.2, 509.2, 511, 516.2, 524.5, 535.4, 548.6, 563.8, 580.6, 598.6, 1321.4, 1339.4,
1356.2, 1371.4, 1384.6, 1395.5, 1403.8, 1409, 1410.8, 1410.8, 1409, 1403.8, 1395.5,
1384.6, 1371.4, 1356.2, 1339.4, 1321.4])
y = np.array([805.4, 805.4, 803.5, 798.3, 790.1,
779.2, 766, 750.8, 734, 716, 364, 346, 329.2, 314, 300.8, 289.9, 281.7, 276.5, 274.6,
274.6, 276.5, 281.7, 289.9, 300.8, 314, 329.2, 346, 364, 716, 734, 750.8, 766, 779.2,
790.1, 798.3, 803.5, 805.4])

fig, ax = plt.subplots(1)
ax.set_aspect('equal')
ax.scatter(x, y, s=40, zorder=3, alpha=0.3)

# compute the distances, ds, between points
dx, dy = x[+1:]-x[:-1],  y[+1:]-y[:-1]
ds = np.array((0, *np.sqrt(dx*dx+dy*dy)))

# compute the total distance from the 1st point, measured on the curve
s = np.cumsum(ds)

# interpolate using 200 point
xinter = np.interp(np.linspace(0,s[-1], 200), s, x)
yinter = np.interp(np.linspace(0,s[-1], 200), s, y)

# plot the interpolated points
ax.scatter(xinter, yinter, s=5, zorder=4)
plt.show()