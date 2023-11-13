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

#NORMALIZATION FACTOR
M = (((0.000000001) / 85 ) / 28 / 1995)

#BINNING
# TAKE HISTOGRAM SIZE DIVIDED BY BINS TO PRODUCE BIN SIZE IN MM
bins = 400

timestr = time.strftime("%Y%m%d-%H%M%S")

beam_c1 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c2 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c3 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c4 = bh.Histogram(bh.axis.Regular((bins),0,2000))

beam_c5 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c6 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c7 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c8 = bh.Histogram(bh.axis.Regular((bins),0,2000))

beam_c9 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c10 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c11 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c12 = bh.Histogram(bh.axis.Regular((bins),0,2000))

beam_c13 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c14 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c15 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c16 = bh.Histogram(bh.axis.Regular((bins),0,2000))

beam_c17 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c18 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c19 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c20 = bh.Histogram(bh.axis.Regular((bins),0,2000))

beam_c21 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c22 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c23 = bh.Histogram(bh.axis.Regular((bins),0,2000))
beam_c24 = bh.Histogram(bh.axis.Regular((bins),0,2000))


name = "1pmt28det"
for rate, hitr, hitp, hitdet, hitpid, hittrid, hitphi in uproot.iterate(["/home/seehrar/projects/def-jmammei/seehrar/outputs/beam_donutremoved_1pmt28det/*.root:T"], ['rate','hit.r','hit.p','hit.det', 'hit.pid', 'hit.trid','hit.ph'], how = tuple, entrysteps = "10MB"):
    #   DETECTOR 28 
    c1  = (hitp < 2) * (hitdet == 28) * ((abs(hitpid) == 11)) * (hitr > 950) * (hitr < 1050) 
    c2  = (hitp > 2) * (hitdet == 28) * ((abs(hitpid) == 11)) * (hitr > 950) * (hitr < 1050)
    c3  = (hitp > 1000) * (hitdet == 28) * ((abs(hitpid) == 11)) * (hitr > 950) * (hitr < 1050)
    c4  = (hitdet == 28) * ((abs(hitpid) == 11)) * (hitr > 950) * (hitr < 1050)
    
    #   DETECTOR 29     
    c5  = (hitp < 2) * (hitdet == 29) * ((abs(hitpid) == 11)) * (hitr > 1100) * (hitr < 1160)
    c6  = (hitp > 2) * (hitdet == 29) * ((abs(hitpid) == 11)) * (hitr > 1100) * (hitr < 1160)
    c7  = (hitp > 1000) * (hitdet == 29) * ((abs(hitpid) == 11)) * (hitr > 1100) * (hitr < 1160)
    c8  = (hitdet == 29) * ((abs(hitpid) == 11)) * (hitr > 1100) * (hitr < 1160)
    
    #   DETECTOR 31   
    c9  = (hitp < 2) * (hitdet == 31) * ((abs(hitpid) == 11)) * (hitr > 1020) * (hitr < 1180)
    c10  = (hitp > 2) * (hitdet == 31) * ((abs(hitpid) == 11)) * (hitr > 1020) * (hitr < 1180)
    c11  = (hitp > 1000) * (hitdet == 31) * ((abs(hitpid) == 11)) * (hitr > 1020) * (hitr < 1180)
    c12  = (hitdet == 31) * ((abs(hitpid) == 11)) * (hitr > 1020) * (hitr < 1180)
    
    #   DETECTOR 32     
    c13  = (hitp < 2) * (hitdet == 31) * ((abs(hitpid) == 11)) * (hitr > 1029) * (hitr < 1189)
    c14  = (hitp > 2) * (hitdet == 31) * ((abs(hitpid) == 11)) * (hitr > 1029) * (hitr < 1189)
    c15  = (hitp > 1000) * (hitdet == 31) * ((abs(hitpid) == 11)) * (hitr > 1029) * (hitr < 1189)
    c16  = (hitdet == 31) * ((abs(hitpid) == 11)) * (hitr > 1029) * (hitr < 1189)
    
    #   DETECTOR 8000   
    c17  = (hitp < 2) * (hitdet == 8000) * ((abs(hitpid) == 11)) * (hitr > 1020) * (hitr < 1180)
    c18  = (hitp > 2) * (hitdet == 8000) * ((abs(hitpid) == 11)) * (hitr > 1020) * (hitr < 1180)
    c19  = (hitp > 1000) * (hitdet == 8000) * ((abs(hitpid) == 11)) * (hitr > 1020) * (hitr < 1180)
    c20  = (hitdet == 8000) * ((abs(hitpid) == 11)) * (hitr > 1020) * (hitr < 1180)
    
    #   DETECTOR 8001     
    c21  = (hitp < 2) * (hitdet == 8001) * ((abs(hitpid) == 11)) * (hitr > 1029) * (hitr < 1189)
    c22  = (hitp > 2) * (hitdet == 8001) * ((abs(hitpid) == 11)) * (hitr > 1029) * (hitr < 1189)
    c23  = (hitp > 1000) * (hitdet == 8001) * ((abs(hitpid) == 11)) * (hitr > 1029) * (hitr < 1189)
    c24  = (hitdet == 8001) * ((abs(hitpid) == 11)) * (hitr > 1029) * (hitr < 1189)


    
    rate_jagged_c1=ak.flatten((ak.ones_like(hitr[c1])) * (rate * M))
    rate_jagged_c2=ak.flatten((ak.ones_like(hitr[c2])) * (rate * M))
    rate_jagged_c3=ak.flatten((ak.ones_like(hitr[c3])) * (rate * M))
    rate_jagged_c4=ak.flatten((ak.ones_like(hitr[c4])) * (rate * M))
    
    rate_jagged_c5=ak.flatten((ak.ones_like(hitr[c5])) * (rate * M))
    rate_jagged_c6=ak.flatten((ak.ones_like(hitr[c6])) * (rate * M))
    rate_jagged_c7=ak.flatten((ak.ones_like(hitr[c7])) * (rate * M))
    rate_jagged_c8=ak.flatten((ak.ones_like(hitr[c8])) * (rate * M))
        
    rate_jagged_c9=ak.flatten((ak.ones_like(hitr[c9])) * (rate * M))
    rate_jagged_c10=ak.flatten((ak.ones_like(hitr[c10])) * (rate * M))
    rate_jagged_c11=ak.flatten((ak.ones_like(hitr[c11])) * (rate * M))
    rate_jagged_c12=ak.flatten((ak.ones_like(hitr[c12])) * (rate * M))
    
    rate_jagged_c13=ak.flatten((ak.ones_like(hitr[c13])) * (rate * M))
    rate_jagged_c14=ak.flatten((ak.ones_like(hitr[c14])) * (rate * M))
    rate_jagged_c15=ak.flatten((ak.ones_like(hitr[c15])) * (rate * M))
    rate_jagged_c16=ak.flatten((ak.ones_like(hitr[c16])) * (rate * M))

    rate_jagged_c17=ak.flatten((ak.ones_like(hitr[c17])) * (rate * M))
    rate_jagged_c18=ak.flatten((ak.ones_like(hitr[c18])) * (rate * M))
    rate_jagged_c19=ak.flatten((ak.ones_like(hitr[c19])) * (rate * M))
    rate_jagged_c20=ak.flatten((ak.ones_like(hitr[c20])) * (rate * M))
    
    rate_jagged_c21=ak.flatten((ak.ones_like(hitr[c21])) * (rate * M))
    rate_jagged_c22=ak.flatten((ak.ones_like(hitr[c22])) * (rate * M))
    rate_jagged_c23=ak.flatten((ak.ones_like(hitr[c23])) * (rate * M))
    rate_jagged_c24=ak.flatten((ak.ones_like(hitr[c24])) * (rate * M))
    
    beam_c1.fill((np.array(ak.flatten(hitr[c1]))), weight = rate_jagged_c1)
    beam_c2.fill((np.array(ak.flatten(hitr[c2]))), weight = rate_jagged_c2)
    beam_c3.fill((np.array(ak.flatten(hitr[c3]))), weight = rate_jagged_c3)
    beam_c4.fill((np.array(ak.flatten(hitr[c4]))), weight = rate_jagged_c4)
    
    beam_c5.fill((np.array(ak.flatten(hitr[c5]))), weight = rate_jagged_c5)
    beam_c6.fill((np.array(ak.flatten(hitr[c6]))), weight = rate_jagged_c6)
    beam_c7.fill((np.array(ak.flatten(hitr[c7]))), weight = rate_jagged_c7)
    beam_c8.fill((np.array(ak.flatten(hitr[c8]))), weight = rate_jagged_c8)
    
    beam_c9.fill((np.array(ak.flatten(hitr[c9]))), weight = rate_jagged_c9)
    beam_c10.fill((np.array(ak.flatten(hitr[c10]))), weight = rate_jagged_c10)
    beam_c11.fill((np.array(ak.flatten(hitr[c11]))), weight = rate_jagged_c11)
    beam_c12.fill((np.array(ak.flatten(hitr[c12]))), weight = rate_jagged_c12)
    
    beam_c13.fill((np.array(ak.flatten(hitr[c13]))), weight = rate_jagged_c13)
    beam_c14.fill((np.array(ak.flatten(hitr[c14]))), weight = rate_jagged_c14)
    beam_c15.fill((np.array(ak.flatten(hitr[c15]))), weight = rate_jagged_c15)
    beam_c16.fill((np.array(ak.flatten(hitr[c16]))), weight = rate_jagged_c16)
    
    beam_c17.fill((np.array(ak.flatten(hitr[c9]))), weight = rate_jagged_c17)
    beam_c18.fill((np.array(ak.flatten(hitr[c10]))), weight = rate_jagged_c18)
    beam_c19.fill((np.array(ak.flatten(hitr[c11]))), weight = rate_jagged_c19)
    beam_c20.fill((np.array(ak.flatten(hitr[c12]))), weight = rate_jagged_c20)
    
    beam_c21.fill((np.array(ak.flatten(hitr[c13]))), weight = rate_jagged_c21)
    beam_c22.fill((np.array(ak.flatten(hitr[c14]))), weight = rate_jagged_c22)
    beam_c23.fill((np.array(ak.flatten(hitr[c15]))), weight = rate_jagged_c23)
    beam_c24.fill((np.array(ak.flatten(hitr[c16]))), weight = rate_jagged_c24)
    
    
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
events_c4 = (np.array(beam_c4)[np.array(beam_c4) != 0])

events_c5 = (np.array(beam_c5)[np.array(beam_c5) != 0])
events_c6 = (np.array(beam_c6)[np.array(beam_c6) != 0])
events_c7 = (np.array(beam_c7)[np.array(beam_c7) != 0])
events_c8 = (np.array(beam_c8)[np.array(beam_c8) != 0])

events_c9 = (np.array(beam_c9)[np.array(beam_c9) != 0])
events_c10 = (np.array(beam_c10)[np.array(beam_c10) != 0])
events_c11 = (np.array(beam_c11)[np.array(beam_c11) != 0])
events_c12 = (np.array(beam_c12)[np.array(beam_c12) != 0])

events_c13 = (np.array(beam_c13)[np.array(beam_c13) != 0])
events_c14 = (np.array(beam_c14)[np.array(beam_c14) != 0])
events_c15 = (np.array(beam_c15)[np.array(beam_c15) != 0])
events_c16 = (np.array(beam_c16)[np.array(beam_c16) != 0])

events_c17 = (np.array(beam_c9)[np.array(beam_c17) != 0])
events_c18 = (np.array(beam_c10)[np.array(beam_c18) != 0])
events_c19 = (np.array(beam_c11)[np.array(beam_c19) != 0])
events_c20 = (np.array(beam_c12)[np.array(beam_c20) != 0])

events_c21 = (np.array(beam_c13)[np.array(beam_c21) != 0])
events_c22 = (np.array(beam_c14)[np.array(beam_c22) != 0])
events_c23 = (np.array(beam_c15)[np.array(beam_c23) != 0])
events_c24 = (np.array(beam_c16)[np.array(beam_c24) != 0])

#NUMBER OF TOTAL EVENTS
count_c1 = np.sum(events_c1 / 8.9e-7)
count_c2 = np.sum(events_c2 / 8.9e-7)
count_c3 = np.sum(events_c3 / 8.9e-7)
count_c4 = np.sum(events_c4 / 8.9e-7)

count_c5 = np.sum(events_c5 / 8.9e-7)
count_c6 = np.sum(events_c6 / 8.9e-7)
count_c7 = np.sum(events_c7 / 8.9e-7)
count_c8 = np.sum(events_c8 / 8.9e-7)

count_c9 = np.sum(events_c9 / 8.9e-7)
count_c10 = np.sum(events_c10 / 8.9e-7)
count_c11 = np.sum(events_c11 / 8.9e-7)
count_c12 = np.sum(events_c12 / 8.9e-7)

count_c13 = np.sum(events_c13 / 8.9e-7)
count_c14 = np.sum(events_c14 / 8.9e-7)
count_c15 = np.sum(events_c15 / 8.9e-7)
count_c16 = np.sum(events_c16 / 8.9e-7)

count_c17 = np.sum(events_c17 / 8.9e-7)
count_c18 = np.sum(events_c18 / 8.9e-7)
count_c19 = np.sum(events_c19 / 8.9e-7)
count_c20 = np.sum(events_c20 / 8.9e-7)

count_c21 = np.sum(events_c21 / 8.9e-7)
count_c22 = np.sum(events_c22 / 8.9e-7)
count_c23 = np.sum(events_c23 / 8.9e-7)
count_c24 = np.sum(events_c24 / 8.9e-7)




df_beam_e = pd.DataFrame()
df_beam_e['Detector'] = ['28','29','31','32']
df_beam_e['(hitp < 2)'] = [("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c1)),(np.sqrt(((8.9e-7) ** 2) * (count_c1))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c5)),(np.sqrt(((8.9e-7) ** 2) * (count_c5))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c9)),(np.sqrt(((8.9e-7) ** 2) * (count_c9))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c13)),(np.sqrt(((8.9e-7) ** 2) * (count_c13))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c17)),(np.sqrt(((8.9e-7) ** 2) * (count_c17))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c21)),(np.sqrt(((8.9e-7) ** 2) * (count_c21)))))]

df_beam_e['(hitp > 2)'] = [("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c2)),(np.sqrt(((8.9e-7) ** 2) * (count_c2))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c6)),(np.sqrt(((8.9e-7) ** 2) * (count_c6))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c10)),(np.sqrt(((8.9e-7) ** 2) * (count_c10))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c14)),(np.sqrt(((8.9e-7) ** 2) * (count_c14))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c18)),(np.sqrt(((8.9e-7) ** 2) * (count_c18))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c22)),(np.sqrt(((8.9e-7) ** 2) * (count_c22)))))]

df_beam_e['(hitp > 1000)'] = [("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c3)),(np.sqrt(((8.9e-7) ** 2) * (count_c3))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c7)),(np.sqrt(((8.9e-7) ** 2) * (count_c7))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c11)),(np.sqrt(((8.9e-7) ** 2) * (count_c11))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c15)),(np.sqrt(((8.9e-7) ** 2) * (count_c15))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c19)),(np.sqrt(((8.9e-7) ** 2) * (count_c19))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c23)),(np.sqrt(((8.9e-7) ** 2) * (count_c23)))))]

df_beam_e['(total)'] = [("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c4)),(np.sqrt(((8.9e-7) ** 2) * (count_c4))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c8)),(np.sqrt(((8.9e-7) ** 2) * (count_c8))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c12)),(np.sqrt(((8.9e-7) ** 2) * (count_c12))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c16)),(np.sqrt(((8.9e-7) ** 2) * (count_c16))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c20)),(np.sqrt(((8.9e-7) ** 2) * (count_c20))))), ("{:1.3f} ± {:1.1e}".format(np.sum(np.array(beam_c24)),(np.sqrt(((8.9e-7) ** 2) * (count_c24)))))]

df_beam_e['UNITS'] = ["GHz/uA/Detector", "GHz/uA/Detector", "GHz/uA/Detector", "GHz/uA/Detector","GHz/uA/Detector","GHz/uA/Detector"]

render_mpl_table(df_beam_e, header_columns=0, col_width=8.5)
plt.savefig("/home/seehrar/projects/def-jmammei/seehrar/analysis/beam_detector_rates_{}_{}.png".format(name, timestr))
