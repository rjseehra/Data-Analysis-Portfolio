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
Detector = 31

# For Detector 31, (hitr > 1020) * (hitr < 1180)
# For Detector 32, (hitr > 1029.059) * (hitr < 1189.057)

#MULTIPLIER
#For everything but pion
M = (((0.000000001) / 85 ) / 14 / 2000)

#BINS
bins = 400

timestr = time.strftime("%Y%m%d-%H%M%S")

beam_c1 = bh.Histogram(bh.axis.Regular((bins),0,2000))
moller_c2 = bh.Histogram(bh.axis.Regular((bins),0,2000))

for rate, hitr, hitp, hitdet, hitpid, hittrid, hitphi in uproot.iterate(["/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/500mill_7cmlucite_showermax_Al_2cmshift_alldet/*.root:T"], ['rate','hit.r','hit.p','hit.det', 'hit.pid', 'hit.trid','hit.ph'], how = tuple, entrysteps = "10MB"):
    c1  = (hitp > 2) * (hitdet == Detector) * ((abs(hitpid) == 11) + (abs(hitpid) == 13) + (abs(hitpid) == 211))
    rate_jagged_c1=ak.flatten((ak.ones_like(hitr[c1])) * (rate * M))
    beam_c1.fill((np.array(ak.flatten(hitr[c1]))), weight = rate_jagged_c1)
    
for rate, hitr, hitp, hitdet, hitpid, hittrid, hitphi in uproot.iterate(["/home/seehrar/projects/rrg-jmammei/elhamm/PionDetectorOptimization/simulation/scratch1/showermaxPlaneAtFrontandBackNoFrameShiftedLucite1cm-WithAluminiumPlate-MoveOutRadially2cm/moller.root:T"], ['rate','hit.r','hit.p','hit.det', 'hit.pid', 'hit.trid','hit.ph'], library="ak", how = tuple):
    c2 = (hitp > 2) * (hitdet == Detector) * ((abs(hitpid) == 11) + (abs(hitpid) == 13) + (abs(hitpid) == 211))
    rate_jagged_c2=ak.flatten((ak.ones_like(hitr[c2])) * (rate/5.95e13))
    moller_c2.fill((np.array(ak.flatten(hitr[c2]))), weight = rate_jagged_c2)
    
plt.figure(1, dpi = 150)
plt.title(" Detector Rates  (Detector {}, Beam)".format(Detector))
plt.ylabel("Particle Rate (GHZ/uA/Detector)")
plt.xlabel("Detector {}".format(Detector))
plt.ylim(0.0000001,10.0)
plt.xlim(0, 2000)
plt.yscale('log')
hep.histplot(beam_c1, color='b', label = "beam_generator")
hep.histplot(moller_c2, color='g', label = "moller generator")
plt.legend(prop={'size': 5}, loc = 'upper right')
plt.savefig("/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/analysis/ratedist_{}_hitp>2_{}.png".format(Detector, timestr))