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

#MULTIPLIER
#For everything but pion
M = (((0.000000001) / 85 ) / 14 / 2000)

#BINS
bins = 400

timestr = time.strftime("%Y%m%d-%H%M%S")

e_p_31 = bh.Histogram(bh.axis.Regular((bins),0,12000))
e_p_32 = bh.Histogram(bh.axis.Regular((bins),0,12000))
muons_31 = bh.Histogram(bh.axis.Regular((bins),0,1000))
muons_32 = bh.Histogram(bh.axis.Regular((bins),0,1000))
pions_31 = bh.Histogram(bh.axis.Regular((bins),0,2000))
pions_32 = bh.Histogram(bh.axis.Regular((bins),0,2000))

for rate, hitr, hitp, hitdet, hitpid, hittrid, hitphi in uproot.iterate(["/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/500mill_7cmlucite_showermax_Al_2cmshift_alldet/*.root:T"], ['rate','hit.r','hit.p','hit.det', 'hit.pid', 'hit.trid','hit.ph'], library="ak", how = tuple):
    c1 = (hitdet == 31) * (abs(hitpid)==11)
    c11 = (hitdet == 32) * (abs(hitpid)==11)
    rate_jagged_c1=ak.flatten((ak.ones_like(hitp[c1])) * (rate * M))
    rate_jagged_c11=ak.flatten((ak.ones_like(hitp[c11])) * (rate * M))
    e_p_31.fill((np.array(ak.flatten(hitp[c1]))), weight = rate_jagged_c1)
    e_p_32.fill((np.array(ak.flatten(hitp[c11]))), weight = rate_jagged_c11)
    
plt.figure(1, dpi = 150)
plt.title(" Energy Spectrum, e + p")
plt.xlabel("Energy (MeV)")
plt.ylabel("Rate (GHz/uA/30 MeV)")
plt.yscale('log')
plt.ylim(0.01, 10000.0)
plt.xlim(0, 12000)
hep.histplot(e_p_31, color='b', label = "e + p @ det 31")
hep.histplot(e_p_32, color='m', alpha = 0.5, label = "e + p @ det 32")
plt.legend()
plt.savefig("/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/analysis/EnergyAnalysis/beam_e_p_{}".format( timestr))

           
for rate, hitr, hitp, hitdet, hitpid, hittrid, hitphi in uproot.iterate(["/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/500mill_7cmlucite_showermax_Al_2cmshift_alldet/*.root:T"], ['rate','hit.r','hit.p','hit.det', 'hit.pid', 'hit.trid','hit.ph'], library="ak", how = tuple):
    c2 = (hitdet == 31) * (abs(hitpid)==13)
    c21 = (hitdet == 32) * (abs(hitpid)==13)
    rate_jagged_c2=ak.flatten((ak.ones_like(hitp[c2])) * (rate * M))
    rate_jagged_c21=ak.flatten((ak.ones_like(hitp[c21])) * (rate * M))
    muons_31.fill((np.array(ak.flatten(hitp[c2]))), weight = rate_jagged_c2)
    muons_32.fill((np.array(ak.flatten(hitp[c21]))), weight = rate_jagged_c21)
    
plt.figure(2, dpi = 150)
plt.title(" Energy Spectrum, muons")
plt.xlabel("Energy (MeV)")
plt.ylabel("Rate (GHz/uA/2.5 MeV)")
plt.yscale('log')
plt.ylim(0.0000001, 0.000001)
plt.xlim(0, 1000)
hep.histplot(muons_31, color='g', label = "muons @ det 31")
hep.histplot(muons_32, color='m', alpha = 0.5, label = "muons @ det 32")
plt.legend()
plt.savefig("/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/analysis/EnergyAnalysis/beam_muons_{}".format(timestr))
           
           
for rate, hitr, hitp, hitdet, hitpid, hittrid, hitphi in uproot.iterate(["/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/500mill_7cmlucite_showermax_Al_2cmshift_alldet/*.root:T"], ['rate','hit.r','hit.p','hit.det', 'hit.pid', 'hit.trid','hit.ph'], library="ak", how = tuple):
    c3 = (hitdet == 31) * (abs(hitpid)==211)
    c31 = (hitdet == 32) * (abs(hitpid)==211)
    rate_jagged_c3=ak.flatten((ak.ones_like(hitp[c3])) * (rate * M))
    rate_jagged_c31=ak.flatten((ak.ones_like(hitp[c31])) * (rate * M))
    pions_31.fill((np.array(ak.flatten(hitp[c3]))), weight = rate_jagged_c3)
    pions_32.fill((np.array(ak.flatten(hitp[c31]))), weight = rate_jagged_c31)
    
plt.figure(3, dpi = 150)
plt.title(" Energy Spectrum, pions")
plt.xlabel("Energy (MeV)")
plt.ylabel("Rate (GHz/uA/5 MeV)")
plt.yscale('log')
plt.ylim(0.0000001, 0.00001)
plt.xlim(0, 2000)
hep.histplot(pions_31, color='y', label = "pions @ det 31")
hep.histplot(pions_32, color='m', alpha = 0.5, label = "pions @ det 32")
plt.legend()
plt.savefig("/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/analysis/EnergyAnalysis/beam_pions_{}".format(timestr))
