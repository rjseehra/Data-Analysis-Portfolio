import sys
#sys.path.append('/home/seehrar/.local/lib/python3.6/site-packages')
#sys.path.append('/home/seehrar/.local/lib/python3.7/site-packages')
#sys.path.append('/home/seehrar/.local/lib/python3.8/site-packages')
import uproot
import numpy as np
import pandas as pd
import matplotlib as mpl 
import matplotlib.pyplot as plt
import awkward1 as ak
import math
from math import pi, sqrt
from matplotlib.collections import LineCollection
from matplotlib import patches
import boost_histogram as bh
from boost_histogram import loc, sum, rebin, underflow, overflow
import time
import tables
import six
import mplhep as hep

#WHAT DETECTOR?
Detector = 8001

#MULTIPLIER
M = (((0.000000001) / 85 ) / 7 / 5000)

#BINS
bins = 400

timestr = time.strftime("%Y%m%d-%H%M%S")

mol_array = bh.Histogram(bh.axis.Regular((bins),24000,26000))
mol_overlay = bh.Histogram(bh.axis.Regular((bins),24000,26000))

name = "moller"
for rate, hitr, hitx, hity, hitz, hitvx, hitvy, hitvz, hitp, hitpx, hitpy, hitpz, hitdet, hitpid, hittrid in uproot.iterate(["moller_*.root"], "T",['rate', 'hit.r','hit.[xyz]', 'hit.vx', 'hit.vy', 'hit.vz', 'hit.p', 'hit.p[xyz]', 'hit.det', 'hit.pid', 'hit.trid'], outputtype = tuple, entrysteps = "10MB"):
    c3 = (np.logical_and(np.logical_and((np.logical_and((hitp > 10), (hitdet == (Detector)))), (hitpid == 11)), (hittrid > 2)))
    c31 = (np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and((np.logical_and((hitp > 10), (hitdet == (Detector)))), (hittrid > 2)), (abs(hitpid) != 2112)), (abs(hitpid) != 22)), (abs(hitpid) != 2212)), (abs(hitpid) != 12)), (abs(hitpid) != 14)))
    plt.figure(6, dpi = 350)
    plt.scatter(hitvz.flatten(),((np.sqrt(hitvx**2+hitvy**2)).flatten()), s = 10, alpha = 0.5, marker = '.', c = 'r', edgecolors = 'black')
plt.title(" sqrt(hitvx**2+hitvy**2) versus hitvz (Detector {}, Moller)".format(Detector), fontsize = 10)
plt.ylabel("Particle Y Position (mm)")
plt.xlabel("Particle Z position (mm)")  
plt.ylim(1000,1500.0)
plt.xlim(24200, 26000)

plt.show()

#plt.yscale('log')
plt.savefig("Tracking/{}_{}_{}".format(name, Detector, timestr), dpi = 150, edgecolor = 'w')


