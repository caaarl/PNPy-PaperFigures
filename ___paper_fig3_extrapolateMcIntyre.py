import matplotlib.pyplot as plt
import numpy as np

# for LaTeX rendering
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True) # doesn't work on Mac, can be uncommented under Linux
rc('text.latex', preamble='\usepackage{sfmath}')

# values from McIntyre et al. 2002
diams = [5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.]
nodeNodeSep = [500, 760, 1000, 1150, 1250, 1350, 1400, 1450, 1500]
noLamella = [80, 100, 110, 120, 130, 135, 140, 145, 150]
nodeDiam = [1.9, 2.4, 2.8, 3.3, 3.7, 4.2, 4.7, 5., 5.5]
MYSADiam = nodeDiam
FLUTLen = [35, 38, 40, 46, 50, 54, 56, 58, 60]
FLUTDiam = [3.4, 4.6, 5.8, 6.9, 8.1, 9.2, 10.4, 11.5, 12.7]
STINLen = [70.5, 111.2, 152.2, 175.2, 190.5, 205.8, 213.5, 221.2, 228.8]
STIDiam = FLUTDiam

diamArray = np.arange(0.1, 20, 0.1)

def cm2inch(value):
    return value/2.54

colors = np.array(((0., 158., 115.), (230., 159., 0.), (86., 180., 233.), (0.,0.,0.)))/255
f, ((ax1, ax2, ax3)) = plt.subplots(1,3, sharex='all',figsize=(8.4,3))

# first plot: diameters
ax1.set_title('diameter')

ax1.scatter(diams, nodeDiam, c=colors[0])
z3 = np.polyfit(diams, nodeDiam, 2)
p3 = np.poly1d(z3)
ax1.plot(diamArray, p3(diamArray), label='Node and MYSA', color=colors[0])

ax1.scatter(diams, FLUTDiam, c=colors[1])
z6 = np.polyfit(diams, FLUTDiam, 2)
p6 = np.poly1d(z6)
ax1.plot(diamArray, p6(diamArray), label='FLUT and STIN', color=colors[1])

ax1.set_ylim(ymin=0)
ax1.set_xlabel('diameter ($\mu$m)')

axbox = ax1.get_position()
ax1.legend(loc = (axbox.x0 - 0.1, axbox.y0 + 0.5), frameon=False)

# second plot: node separation
ax2.scatter(diams, nodeNodeSep, c=colors[0])
ax2.set_title('node separation')
z1 = np.polyfit(diams, nodeNodeSep, 1)
p1 = np.poly1d(z1)
ax2.plot(diamArray, p1(diamArray), color=colors[0])
ax2.set_xlabel('diameter ($\mu$m)')

# third plot: lamella count
ax3.scatter(diams, noLamella, c=colors[0])
ax3.set_title('number of lamella')
z2 = np.polyfit(diams, noLamella, 1)
p2 = np.poly1d(z2)
ax3.plot(diamArray, p2(diamArray), color=colors[0])
ax3.set_xlim([0,20])
ax3.set_xlabel('diameter ($\mu$m)')

# remove unneccessary axis lines
for ax in (ax1, ax2, ax3):
    for loc, spine in ax.spines.items():
        if loc in ['right', 'top']:
            spine.set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

# label panels
panelLabels = ['A', 'B', 'C']
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_weight('heavy')

for axInd, ax in enumerate((ax1, ax2, ax3)):
    plt.text(0.1, 0.9, panelLabels[axInd],
         horizontalalignment='center',
         verticalalignment='center',
         fontsize=12,
         fontweight='heavy',
         fontproperties=font,
         transform = ax.transAxes)

plt.tight_layout()

import os
if not os.path.exists('figures'):
    os.makedirs('figures')
plt.savefig(os.path.join('figures', 'fig3_extrapolateMcIntyre.eps'), format='eps', dpi=300)

plt.show()

