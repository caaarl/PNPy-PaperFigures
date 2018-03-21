import cPickle as pickle
import os
import numpy as np
import matplotlib.pyplot as plt

colors = np.array(((0.,0.,0.), (230., 159., 0.), (86., 180., 233.), (0., 158., 115.)))/255

# load SFAP dictionary
saveDict = pickle.load(open(os.path.join('SFAPs', 'SFAPsCuffWidth3.dict'), "rb" ))

unmyelinatedCuffWidths = saveDict['cuffWidthsUnmyelinated']
myelinatedCuffWidths = saveDict['cuffWidthsMyelinated']
unmyelinatedSFAPs = saveDict['unmyelinatedSFAPs']
myelinatedSFAPs = saveDict['myelinatedSFAPs']
diametersUnmyel = saveDict['unmyelinatedDiameters']
diametersMyel = saveDict['myelinatedDiameters']

amplitudesMyel = np.zeros((len(diametersMyel), len(myelinatedCuffWidths)))
amplitudesUnmyel = np.zeros((len(diametersUnmyel), len(unmyelinatedCuffWidths)))

for diamInd in range(len(diametersUnmyel)):
    for cuffWidthInd in range(len(unmyelinatedCuffWidths)):
        amplitudesUnmyel[diamInd, cuffWidthInd] = np.max(unmyelinatedSFAPs[diamInd][cuffWidthInd]) - \
                                                  np.min(unmyelinatedSFAPs[diamInd][cuffWidthInd])

for diamInd in range(len(diametersMyel)):
    for cuffWidthInd in range(len(myelinatedCuffWidths)):
        amplitudesMyel[diamInd, cuffWidthInd] = np.max(myelinatedSFAPs[diamInd][cuffWidthInd]) - \
                                                np.min(myelinatedSFAPs[diamInd][cuffWidthInd])

f, axarr = plt.subplots(1,2)

axarr[0].imshow(np.multiply(amplitudesUnmyel.T,1000), interpolation='bilinear')
c0 = axarr[0].contour(np.multiply(amplitudesUnmyel.T,1000), colors='w')
plt.clabel(c0, inline=1, fontsize=10)

axarr[1].imshow(np.multiply(amplitudesMyel.T,1000), interpolation='bilinear')
c1 = axarr[1].contour(np.multiply(amplitudesMyel.T,1000), colors='w')
plt.clabel(c1, inline=1, fontsize=10)

diameterArrays = [diametersUnmyel, diametersMyel]
cuffWidths = [unmyelinatedCuffWidths, myelinatedCuffWidths]
for i in range(2):
    axarr[i].set_yticks(range(len(cuffWidths[i])))
    yTicks = []
    for cuffWidth in cuffWidths[i]:
        yTicks.append('%1.2f' % np.multiply(cuffWidth,1000))
    axarr[i].set_yticklabels(yTicks)

    axarr[i].set_ylabel('cuff width (mm)')

    axarr[i].set_xticks(range(len(diameterArrays[i])))
    axarr[i].set_xticklabels(diameterArrays[i], rotation=45)
    axarr[i].set_xlabel('diameter (um)')

if not os.path.exists('figures'):
    os.makedirs('figures')
    
plt.savefig(os.path.join('figures', 'fig9_cuffWidthAmplitude.eps'),
        format='eps', dpi=300)

plt.show()