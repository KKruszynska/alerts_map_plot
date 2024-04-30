import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
from matplotlib.legend_handler import HandlerTuple

import sys

import astropy.coordinates as coord
import astropy.units as u

from imageio import imread

datafile = "Gaia_s_sky_in_colour_mod3.png"
img = imread(datafile)

#gsa_list = "/home/kasia/Documents/PhD/Gaia/gaia_alerts/Gaia_list.dat"
#gsa_list = "/home/kasia/Documents/PhD/Gaia/gaia_alerts/Gaia_list1620.dat"
# DR3_listA1 = "/home/kasia/Documents/PhD/Gaia/cu7/SOS/MicrolensingAstro/coordinate_list_A1.csv"
# DR3_listB1 = "/home/kasia/Documents/PhD/Gaia/cu7/SOS/MicrolensingAstro/coordinate_list_B1.csv"

DR3_listA1 = "/home/kasia/Documents/PhD/Gaia/cu7/ForPaper/DR3_vari_table_final_merged_corr_lvl2_names.csv"

#gaia_name = np.genfromtxt(gsa_list, usecols=0,dtype=str, unpack=True)
#gaia_l, gaia_b = np.genfromtxt(gsa_list, usecols=(1,2), dtype=np.float, unpack=True)

listA1 = pd.read_csv(DR3_listA1, header=0)
# listB1 = pd.read_csv(DR3_listB1, header=0)

galA1 = coord.SkyCoord(ra=listA1["ra_deg"]*u.degree, dec=listA1["dec_deg"]*u.degree)
# galB1 = coord.SkyCoord(ra=listB1["ra_deg"]*u.degree, dec=listB1["dec_deg"]*u.degree)

# for i in range(len(galA1)):
	# if(galA1[i].galactic.b.degree>15.):
		# print(listA1["sourceid_1"].values[i])
		
# Recalculating the coordinates
lg = coord.Angle(galA1.galactic.l.degree*u.degree)
lg = lg.wrap_at(180*u.degree)
lg = (-1)*lg
bg = coord.Angle(galA1.galactic.b.degree*u.degree)

# ll = coord.Angle(galB1.galactic.l.degree*u.degree)
# ll = ll.wrap_at(180*u.degree)
# ll = (-1)*ll
# bl = coord.Angle(galB1.galactic.b.degree*u.degree)

# lgsa = coord.Angle(gaia_l*u.degree)
# lgsa = lgsa.wrap_at(180*u.degree)
# lgsa = (-1)*lgsa
# bgsa = coord.Angle(gaia_b*u.degree)

fig = plt.figure(figsize=(15,7))
ax0 = fig.add_subplot(111, label="Cartesian background")
ax0.imshow(img, origin="upper", zorder=0)
ax0.axis("off")
ax = fig.add_subplot(111, projection="aitoff", label="Galaxy map")
ax.set_facecolor("None")
ax.tick_params(axis='x', colors='white', direction='in')
#ax.tick_params(axis='y', colors='white', direction='in')
ax.tick_params(axis='y', colors='black', direction='in')

gaiaA1 = ax.scatter(lg.radian,bg.radian, c='yellow', edgecolor='black', marker='o',s=80, zorder=0)
# gaiaB1 = ax.scatter(ll.radian,bl.radian, c='yellow', edgecolor='black', marker='^',s=80, zorder=0)
#gaiaGSA = ax.scatter(lgsa.radian,bgsa.radian, c='yellow', edgecolor='black', marker='o',s=80, zorder=0)

ax.set_xticklabels(['150$^\circ$','120$^\circ$','90$^\circ$','60$^\circ$','30$^\circ$','0$^\circ$','330$^\circ$','300$^\circ$','270$^\circ$','240$^\circ$','210$^\circ$'])
ax.grid(True)

print()
# label1 = 'Gaia DR3 ('+str(len(listA1["ra_deg"])+len(listB1["ra_deg"])-33)+')'
label1 = 'Gaia DR3 ('+str(len(listA1["ra_deg"]))+')'
#label2 = 'Gaia Science Alerts ('+str(len(gaia_name))+')'

# plt.legend([(gaiaA1, gaiaB1)], [label1], bbox_to_anchor=(0.25, 0.05), fontsize=16, handler_map={tuple: HandlerTuple(ndivide=None)})
plt.legend([(gaiaA1)], [label1], bbox_to_anchor=(0.25, 0.05), fontsize=16, handler_map={tuple: HandlerTuple(ndivide=None)})
#plt.legend([gaiaGSA], [label2], bbox_to_anchor=(0.25, 0.05), fontsize=16)
#plt.legend([gaia, con, ogle, asassn], [label1,label2,label3,label4], bbox_to_anchor=(0.25, 0.05), fontsize=14)
#plt.legend([gaiaA1, gaiaGSA], [label1, label2], bbox_to_anchor=(0.25, 0.05), fontsize=14)
plt.savefig("Gaia_alerts_dr3_2024_jan.png", transparent=True)
plt.show()
