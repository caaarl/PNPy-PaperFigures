import cPickle as pickle
import os
import numpy as np
import matplotlib.pyplot as plt

colors = np.array(((0.,0.,0.), (230., 159., 0.), (86., 180., 233.), (0., 158., 115.)))/255

# 'SFAPsPowleyMyelAsRecordingsIdealizedCuffHigherDiams.dict'
saveDict = pickle.load(open(os.path.join('SFAPs', 'SFAPsForAmplitudePlot.dict'), 'rb')) #'SFAPsPowleyMyelAsRecordingsIdealizedCuff2.dict'), "rb" )) # originalMyelDiam #/Volumes/SANDISK/PyPN/

stringsDiam = ['unmyelinatedDiameters', 'myelinatedDiameters']
stringsSFAPHomo = ['unmyelinatedSFAPsHomo', 'myelinatedSFAPsHomo']
stringsSFAPFEM = ['unmyelinatedSFAPsFEM', 'myelinatedSFAPsFEM']
stringsSFAPIdeal = ['unmyelinatedSFAPIdeal', 'myelinatedSFAPIdeal']

fieldStrings = ['homogeneous', 'radially inhomogeneous', 'cuff']
axonTypeStrings = ['unmyelinated', 'myelinated']
lineTypes = ['-', ':']

# medium types
for fieldTypeInd in [2,1,0]: # fieldTypes:

    # axon types
    for typeInd in [0,1]:

        diameters = np.array(saveDict[stringsDiam[typeInd]])

        if fieldTypeInd == 0:
            SFAP = np.transpose(np.array(saveDict[stringsSFAPHomo[typeInd]]))
        elif fieldTypeInd == 1:
            SFAP = np.transpose(np.array(saveDict[stringsSFAPFEM[typeInd]]))
        else:
            SFAP = np.transpose(np.array(saveDict[stringsSFAPIdeal[typeInd]]))

        amplitudes = np.max(SFAP, 0) - np.min(SFAP, 0)

        plt.semilogy(diameters, amplitudes, lineTypes[typeInd], color=colors[fieldTypeInd+1], linewidth=2,
                     label=fieldStrings[fieldTypeInd]+' '+axonTypeStrings[typeInd])


plt.legend(loc='best', frameon=False)
plt.xlabel('diameter ($\mu$m)')
plt.ylabel('amplitude ($\mu$V)')

if not os.path.exists('figures'):
    os.makedirs('figures')
    
plt.savefig(os.path.join('figures', 'SFAPAmps.eps'),
        format='eps', dpi=300)

plt.show()