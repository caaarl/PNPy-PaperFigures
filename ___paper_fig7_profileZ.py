import PyPNS
import numpy as np
import matplotlib.pyplot as plt

# for LaTeX rendering
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)
rc('text.latex', preamble='\usepackage{sfmath}')

# set length of bundle and number of axons
lengthOfBundle = 400000
segmentLengthAxon = 30
bundleGuide = PyPNS.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

extracellulars = []
extracellulars.append(PyPNS.Extracellular.homogeneous(sigma=1.2))
extracellulars.append(PyPNS.Extracellular.precomputedFEM(bundleGuide))
extracellulars.append(PyPNS.Extracellular.analytic(bundleGuide))

xPositions = np.arange(-15000, 15000, 10)
sourceCurrents = np.ones(1)

colors = np.array(((0.,0.,0.), (230., 159., 0.), (86., 180., 233.), (0., 158., 115.)))/255
lineStyles = ['-', '--', '-.']

legends = ['hom.', 'rad. inhom.', 'cuff']
for extraInd, extracellular in enumerate(extracellulars):

    receiverPositions = np.hstack((xPositions[:, np.newaxis], np.zeros((len(xPositions), 1)), np.ones((len(xPositions), 1)) * 235))

    for xPInd, xP in enumerate([0, 180]):

        sourcePositions = np.array([0, 0, xP])

        v = extracellular.calculate_extracellular_potential(np.transpose(sourcePositions[:, np.newaxis]), sourceCurrents[:, np.newaxis], receiverPositions)
        plt.plot(xPositions, v / np.max(v), lineStyles[xPInd], label=legends[extraInd] + ', xP=' + str(xP) + ' $\mu$m', color=colors[extraInd + 1], linewidth=2.0)

plt.ylabel('normalised amplitude')
plt.xlabel('longitudinal distance ($\mu$m)')
plt.legend(loc='best', frameon=False)

import os
if not os.path.exists('figures'):
    os.makedirs('figures')
    
xlimits = ([-15000,15000], [-500,500])
figureNames = ['fig7_voltageProfileFull.eps', 'fig7_voltageProfileZoomed.eps']
for limInd in range(2):
    plt.xlim(xlimits[limInd])

    plt.savefig(os.path.join('figures', figureNames[limInd]),
            format='eps', dpi=300)

plt.show()