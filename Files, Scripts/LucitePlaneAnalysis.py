#!/usr/bin/env python
import sys
#sys.path.append('/home/seehrar/.local/lib/python3.6/site-packages')
#sys.path.append('/home/seehrar/.local/lib/python3.7/site-packages')
#sys.path.append('/home/seehrar/.local/lib/python3.8/site-packages')
import uproot
import numpy as np
import pandas as pd
import matplotlib as mpl 
import matplotlib.pyplot as plt
import awkward as ak
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

timestr = time.strftime("%Y%m%d-%H%M%S")

name = "7cmlucite_S+Al"
for hitr, hitx, hity, hitz, hitp, hitdet, hitpid, hittrid in uproot.iterate(["/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/500mill_7cmlucite_showermax_Al_1cmshift/*.root"],['hit.r','hit.x','hit.y','hit.z','hit.p','hit.det','hit.pid','hit.trid'], how = tuple, entrysteps = "10MB"):
    c1 = (hitdet == 8000) * (hitp < 10e-6) #npes
    c2 = (hitdet == 8001) * ((abs(hitpid) == 11) + (abs(hitpid) == 13) + (abs(hitpid) == 211)) 
    plt.figure(1, dpi=150)
    plt.scatter(ak.flatten(hitx[c1]), ak.flatten(hity[c1]), s = 10, alpha = 0.5, marker = '.', c = 'r', edgecolors = 'red')
    plt.scatter(ak.flatten(hitx[c2]), ak.flatten(hity[c2]), s = 10, alpha = 0.5, marker = '.', c = 'b', edgecolors = 'blue')
plt.title("Lucite XY-Plane Analysis", fontsize = 10)
plt.ylabel("Particle Y Position (mm)")
plt.xlabel("Particle X position (mm)")
plt.ylim(-1300,-900)
plt.xlim(0, 510)
plt.legend(["Detector 8000", "Detector 8001"], loc ="upper right")
plt.show()
#plt.yscale('log')
plt.savefig("LuciteXYAnalysis/{}_{}.png".format(name,timestr), dpi = 150, edgecolor = 'w')
