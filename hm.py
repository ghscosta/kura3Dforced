import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
from mpl_toolkits import mplot3d
from numpy import linalg as LA
import pandas as pd
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import matplotlib 
import os
from glob import glob as g
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import math as math

dados = pd.read_csv("fig5New.dat",sep='\s+')
#dados = pd.read_csv("bif-1-F-ome-n.csv")

#dados = dados.query('sigma <= 1.0')
#dados = dados.query('force <= 1.0')

f,(ax1) = plt.subplots(1,1, figsize=(9,9))
cbar_ax = f.add_axes([.15, -0.03, .7, .04])
cbarticks = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
matplotlib.rcParams["mathtext.fontset"] = 'cm'
coluna = 'rmodav'
pivotados = dados.pivot(index='force',columns='sigma',values=coluna)

sns.heatmap(pivotados, square = False, ax=ax1,
            cmap = "viridis",zorder=1,xticklabels=20,yticklabels=80,cbar_ax = cbar_ax,
            cbar_kws={'orientation':"horizontal","ticks":cbarticks})

#ax1.set_yticklabels(yticklabels,rotation=0)
ax1.tick_params(labelsize=20)
ax1.set_xlabel(r"$K_1$",fontsize=30)
ax1.set_ylabel(r"$K_2$",fontsize=30)
ax1.annotate("(c)", xy=(0.1,0.8), xytext=(60,40),fontsize = 30,color='black')


cax = ax1.figure.axes[-1]
cax.tick_params(labelsize=20)
cax.set_frame_on(True)
ax1.invert_yaxis()
f.savefig("Fig5"+coluna+".png",bbox_extra_artists=(cbar_ax,),bbox_inches='tight')
f.savefig("Fig5"+coluna+".pdf",bbox_extra_artists=(cbar_ax,),bbox_inches='tight')