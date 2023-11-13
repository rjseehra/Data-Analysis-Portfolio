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
Detector = 32

# For Detector 31, (hitr > 1020) * (hitr < 1180)
# For Detector 32, (hitr > 1029.059) * (hitr < 1189.057)

#MULTIPLIER
#For everything but pion
M = (((0.000000001) / 85 ) / 14 / 2000)

#BINS
bins = 400

timestr = time.strftime("%Y%m%d-%H%M%S")

beam_c1 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c2 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c3 = bh.Histogram(bh.axis.Regular((bins),0,2000))

name = "beam"
for rate, hitr, hitp, hitdet, hitpid, hittrid, hitphi in uproot.iterate(["/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/500mill_7cmlucite_showermax_Al_2cmshift_alldet/*.root:T"], ['rate','hit.r','hit.p','hit.det', 'hit.pid', 'hit.trid','hit.ph'], how = tuple, entrysteps = "10MB"):
    c1  = (hitp < 2) * (hitdet == Detector) * ((abs(hitpid) == 11) + (abs(hitpid) == 13) + (abs(hitpid) == 211)) * (hitr > 1029.059) * (hitr < 1189.059)
    c2  = (hitp > 2) * (hitdet == Detector) * ((abs(hitpid) == 11) + (abs(hitpid) == 13) + (abs(hitpid) == 211)) * (hitr > 1029.059) * (hitr < 1189.059)
    c3  = (hitdet == Detector) * ((abs(hitpid) == 11) + (abs(hitpid) == 13) + (abs(hitpid) == 211)) * (hitr > 1029.059) * (hitr < 1189.059)
    
    rate_jagged_c1=ak.flatten((ak.ones_like(hitr[c1])) * (rate * M))
    rate_jagged_c2=ak.flatten((ak.ones_like(hitr[c2])) * (rate * M))
    rate_jagged_c3=ak.flatten((ak.ones_like(hitr[c3])) * (rate * M))
    
    beam_c1.fill((np.array(ak.flatten(hitr[c1]))), weight = rate_jagged_c1)
    beam_c2.fill((np.array(ak.flatten(hitr[c2]))), weight = rate_jagged_c2)
    beam_c3.fill((np.array(ak.flatten(hitr[c3]))), weight = rate_jagged_c3)
    
def render_mpl_table(data, col_width=1.0, row_height=1.625, font_size=27,
                     header_color='b', row_colors=['w', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 4])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')

    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in  six.iteritems(mpl_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    return ax

#CREATING EVENT ARRAY
events_c1 = (np.array(beam_c1)[np.array(beam_c1) != 0])
events_c2 = (np.array(beam_c2)[np.array(beam_c2) != 0])
events_c3 = (np.array(beam_c3)[np.array(beam_c3) != 0])
#NUMBER OF TOTAL EVENTS
count_c1 = np.sum(events_c1 / 8.9e-7)
count_c2 = np.sum(events_c2 / 8.9e-7)
count_c3 = np.sum(events_c3 / 8.9e-7)

df_beam_c1 = pd.DataFrame()
df_beam_c1['Gen. ({})'.format(Detector)] = ['{}'.format(name)]
df_beam_c1['Rate (hitp < 2)'] = [("{:1.3e} ± {:1.1e}".format(np.sum(np.array(beam_c1)),(np.sqrt(((8.9e-7) ** 2) * (count_c1)))))]
df_beam_c1['UNITS'] = ["GHz/uA/Detector"]

df_beam_c2 = pd.DataFrame()
df_beam_c2['Gen. ({})'.format(Detector)] = ['{}'.format(name)]
df_beam_c2['Rate (hitp > 2)'] = [("{:1.3e} ± {:1.1e}".format(np.sum(np.array(beam_c2)),(np.sqrt(((8.9e-7) ** 2) * (count_c2)))))]
df_beam_c2['UNITS'] = ["GHz/uA/Detector"]

df_beam_c3 = pd.DataFrame()
df_beam_c3['Gen. ({})'.format(Detector)] = ['{}'.format(name)]
df_beam_c3['Rate (total)'] = [("{:1.3e} ± {:1.1e}".format(np.sum(np.array(beam_c3)),(np.sqrt(((8.9e-7) ** 2) * (count_c3)))))]
df_beam_c3['UNITS'] = ["GHz/uA/Detector"]

render_mpl_table(df_beam_c1, header_columns=0, col_width=6.5)
plt.savefig("/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/analysis/beam_{}_hitp_less_than_2_{}.png".format(Detector, timestr))
render_mpl_table(df_beam_c2, header_columns=0, col_width=6.5)
plt.savefig("/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/analysis/beam_{}_hitp_greater_than_2_{}.png".format(Detector, timestr))
render_mpl_table(df_beam_c3, header_columns=0, col_width=6.5)
plt.savefig("/home/seehrar/projects/rrg-jmammei/seehrar/simulation_analysis/cedar_2.0/analysis/beam_{}_total_{}.png".format(Detector, timestr))
