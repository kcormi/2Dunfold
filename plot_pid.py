import numpy as np
from argparse import ArgumentParser
import os
import matplotlib as mpl
mpl.use('Agg')
import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.ticker import MaxNLocator
hep.style.use("CMS")
import ast

particle = ["$e^-$","$\pi^+$","$K^+$","p","$\Sigma$","$\Xi^-$","$\Omega^-$"]
N_EPOS = np.array([233,179503,33854,17498,2634,819,91])
N_CP5 = np.array([376,169961,30126,15017,2042,524,13])
N_CP1 = np.array([653,258547,44635,24318,3123,933,33])
e_EPOS = np.sqrt(N_EPOS)
e_CP5 = np.sqrt(N_CP5)
e_CP1 = np.sqrt(N_CP1)
ratio_EPOS = N_EPOS/sum(N_EPOS)
e_ratio_EPOS = e_EPOS/sum(N_EPOS)
ratio_CP5 = N_CP5/sum(N_CP5)
e_ratio_CP5 = e_CP5/sum(N_CP5)
ratio_CP1 = N_CP1/sum(N_CP1)
e_ratio_CP1 = e_CP1/sum(N_CP1)
reldiff_CP5_EPOS = (ratio_CP5-ratio_EPOS)/ratio_EPOS
e_reldiff_CP5_EPOS = np.sqrt(np.square(e_ratio_CP5/ratio_EPOS)+np.square(ratio_CP5*e_ratio_EPOS/np.square(ratio_EPOS)))
reldiff_CP1_EPOS = (ratio_CP1-ratio_EPOS)/ratio_EPOS
e_reldiff_CP1_EPOS = np.sqrt(np.square(e_ratio_CP1/ratio_EPOS)+np.square(ratio_CP1*e_ratio_EPOS/np.square(ratio_EPOS)))

color = ['r' if d>0 else 'b' for d in reldiff_CP5_EPOS]
plt.bar(range(len(particle)),reldiff_CP5_EPOS,yerr=e_reldiff_CP5_EPOS,color=color)
plt.xticks(range(len(particle)),particle)
plt.title("relative difference of Nch ratio CP5 versus EPOS",{'fontsize':20})
plt.savefig("piddiff_CP5_EPOS.pdf")
plt.close()

color = ['r' if d>0 else 'b' for d in reldiff_CP1_EPOS]
plt.bar(range(len(particle)),reldiff_CP1_EPOS,yerr=e_reldiff_CP1_EPOS,color=color)
plt.xticks(range(len(particle)),particle)
plt.title("relative difference of Nch ratio CP1 versus EPOS",{'fontsize':20})
plt.savefig("piddiff_CP1_EPOS.pdf")
plt.close()

plt.bar(np.arange(len(particle))-0.25,ratio_EPOS,width=0.25,yerr=e_ratio_EPOS,label="EPOS")
plt.bar(np.arange(len(particle)),ratio_CP5,width=0.25,yerr=e_ratio_CP5,label="CP5")
plt.bar(np.arange(len(particle))+0.25,ratio_CP1,width=0.25,yerr=e_ratio_CP1,label="CP1")
plt.xticks(range(len(particle)),particle)
plt.title("Nch ratios",{'fontsize':20})
plt.legend(loc='upper right')
plt.yscale("log")
plt.savefig("pid_CP5_CP1_EPOS.pdf")
plt.close()
