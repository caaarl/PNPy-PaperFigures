import cPickle as pickle
import os
import numpy as np
import matplotlib.pyplot as plt

colors = np.array(((0.,0.,0.), (230., 159., 0.), (86., 180., 233.), (0., 158., 115.)))/255

# 'SFAPsPowleyMyelAsRecordingsIdealizedCuffHigherDiams.dict'
saveDict = pickle.load(open(os.path.join('iMems', 'iMemIntVsDiam3.dict'), 'rb')) #'SFAPsPowleyMyelAsRecordingsIdealizedCuff2.dict'), "rb" )) # originalMyelDiam #/Volumes/SANDISK/PyPN/

stringsDiam = ['unmyelinatedDiameters', 'myelinatedDiameters']
stringsIMem = ['unmyelinatedIMem', 'myelinatedIMem']

axonTypeStrings = ['unmyelinated', 'myelinated']
lineTypes = ['-', ':']

# axon types
for typeInd in [0,1]:

    diameters = np.array(saveDict[stringsDiam[typeInd]])

    iMemInts = saveDict[stringsIMem[typeInd]]

    plt.semilogy(diameters, iMemInts, color=colors[typeInd+1], linewidth=2, # , lineTypes[typeInd]
                 label=axonTypeStrings[typeInd])


plt.legend(loc='best', frameon=False)
plt.xlabel('diameter ($\mu$m)')
plt.ylabel('integrated current (nA ms / $\mu$m)')

if not os.path.exists('figures'):
    os.makedirs('figures')
    
plt.savefig(os.path.join('figures', 'iMemInt.eps'),
        format='eps', dpi=300)

plt.show()