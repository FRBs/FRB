#!/usr/bin/env python
import os, sys
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

class Localisation:
    def __init__(self, rastring, decstring, uncertaintymaj, uncertaintymin, uncertaintypadeg):
        self.pos = SkyCoord(rastring + " " + decstring)
        self.uncertaintymaj = uncertaintymaj
        self.uncertaintymin = uncertaintymin
        self.uncertaintypadeg = uncertaintypadeg

    def chisq(self, position):
        offset = position.separation(self.pos).arcsecond
        pa = position.position_angle(self.pos).degree
        deltaparad = (pa - self.uncertaintypadeg)*np.pi/180
        uncertainty = np.sqrt((self.uncertaintymaj*np.cos(deltaparad))**2 + (self.uncertaintymin*np.sin(deltaparad))**2)
        #print(self.pos, position, offset, uncertainty)
        return (offset / uncertainty)**2

# The map centre is the VLBI localisation
mapcentre = SkyCoord("05h08m03.5077s +26d03m38.504s")
constraints = []
constraints.append(Localisation("05h08m03.475s", "+26d03m38.44s", 1.25, 0.67, 140))
constraints.append(Localisation("05h08m03.665s", "+26d03m39.50s", 1.55, 0.91, 42))
constraints.append(Localisation("05h08m03.540s", "+26d03m39.05s", 1.14, 0.45, 0))

gridsize = 50
gridspacing = 0.1 #arcsec
chisqgrid = np.zeros(gridsize*gridsize).reshape(gridsize, gridsize)
x = np.zeros(gridsize*gridsize).reshape(gridsize, gridsize)
y = np.zeros(gridsize*gridsize).reshape(gridsize, gridsize)

for xidx in range(gridsize):
    xoffset = (xidx - gridsize//2)*gridspacing
    print("Working on row",xidx+1,"/",gridsize)
    for yidx in range(gridsize):
        yoffset = (yidx - gridsize//2)*gridspacing
        x[xidx, yidx] = xoffset
        y[xidx, yidx] = yoffset
        sc = SkyCoord(mapcentre.ra + xoffset*u.arcsec, mapcentre.dec + yoffset*u.arcsec)
        for c in constraints:
            chisqgrid[xidx][yidx] += c.chisq(sc)
pgrid = np.exp(-0.5*np.square(chisqgrid))
pgrid /= np.sum(pgrid)
ex = [x[0,0], x[gridsize-1,gridsize-1], y[0,0], y[gridsize-1,gridsize-1]]

print("Minimum chi squared is ", np.min(chisqgrid))
m = np.min(chisqgrid)

cmap = plt.get_cmap('PiYG')
#levels = MaxNLocator(nbins=3).tick_values(1, 10)
levels = [m, m+2.3, m+6.17, m+11.8]
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
cf = plt.contourf(x, y, chisqgrid, levels=levels, cmap=cmap)
plt.savefig("chisqplot.png")

plt.clf()
n = 1000
t = np.linspace(0, pgrid.max(), n)
integral = ((pgrid >= t[:, None, None]) * pgrid).sum(axis=(1,2))
from scipy import interpolate
f = interpolate.interp1d(integral, t)
t_contours = f(np.array([0.95, 0.68]))

plt.imshow(pgrid, origin='lower', extent=ex, cmap="gray")
plt.contour(pgrid, t_contours, extent=ex)
plt.savefig("probplot.png")

