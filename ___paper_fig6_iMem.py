import PyPNS
import matplotlib.pyplot as plt
import numpy as np

# axon definitions
myelinatedParameters = {'fiberD': 3} # , 'gkbar_axnode': gkbar} # , 'temperature': temperature}
unmyelinatedParameters = {'fiberD': 3} # , 'temperature': temperature}

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 10., # 100. #50,  # Pulse amplitude (mA)
                           'frequency': 20.,  # Frequency of the pulse (kHz)
                           'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                           'stimDur': 0.05,  # Stimulus duration (ms)
                           'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                           'delay': 0.,  # ms
                           }
intraParameters = {'stimulusSignal': PyPNS.signalGeneration.rectangular(**rectangularSignalParams)}

# plot configuration; colorblind compatible colors
colors = np.array(((0., 158., 115.), (230., 159., 0.), (86., 180., 233.), (0.,0.,0.)))/255
(f,axarr) = plt.subplots(1,2,sharey='all',figsize=(8.4, 4))

# loop over both axon types; 0: unmyelinated, 1: myelinated
for typeInd in [0,1]:

    # set all properties of the bundle
    bundleParameters = {'radius': 200,  #um Radius of the bundle
                        'length': 5000,  # um Axon length

                        'numberOfAxons': 1,  # Number of axons in the bundle
                        'pMyel': typeInd,  # Percentage of myelinated fiber type A
                        'pUnmyel': 1 - typeInd,  # Percentage of unmyelinated fiber type C
                        'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                        'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                        'tStop': 100, # ms
                        'timeRes': 0.0025, # ms

                        'saveI': True, # save current (False by default)
                        'saveV': False, # don't save voltage
                        }

    # create the bundle with all properties of axons and recording setup
    bundle = PyPNS.Bundle(**bundleParameters)

    # spiking through a single electrical stimulation
    bundle.add_excitation_mechanism(PyPNS.StimIntra(**intraParameters))

    # run the simulation
    bundle.simulate()

    # load current
    t, imem = bundle.get_imem_from_file_axonwise(0)

    if typeInd == 1:

        nodeIndices = np.where(np.array(bundle.axons[0].segmentTypes) == 'n')[0]
        nodeIndex = nodeIndices[int(len(nodeIndices) * 0.8)]

        labelStringsSections = ['node', 'MYSA', 'FLUT']
        for ii in range(3):
            axarr[typeInd].plot(t, imem[nodeIndex + ii, :], label=labelStringsSections[ii], color=colors[ii+1,:], linewidth=2)
        axarr[typeInd].plot(t, np.sum(imem[nodeIndex:nodeIndex + 3, :], axis=0), label='sum of all',
                            linewidth=2, color=colors[0,:])

        axarr[typeInd].set_title('Myelinated')
        axarr[typeInd].legend(loc='best', frameon=False)

        axarr[typeInd].set_xlim((0.2, 1.3))
    else:
        nodeIndex = 10
        axarr[typeInd].plot(t, imem[nodeIndex, :], color=colors[0,:], linewidth=2)
        axarr[typeInd].set_title('Unmyelinated')
        axarr[typeInd].set_ylabel('$I_m$ (nA)')

        axarr[typeInd].set_ylim((-1.6, 0.75))
        axarr[typeInd].set_xlim((0, 1.5))

    axarr[typeInd].set_xlabel('time (ms)')

# remove unneccessary axis lines
for ax in axarr:
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

for axInd, ax in enumerate(axarr):
    plt.text(0.1, 0.9, panelLabels[axInd],
         horizontalalignment='center',
         verticalalignment='center',
         fontsize=12,
         fontweight='heavy',
         fontproperties=font,
         transform = ax.transAxes)

import os
if not os.path.exists('figures'):
    os.makedirs('figures')
plt.savefig(os.path.join('figures', 'fig6_imemTimecourse.eps'), format='eps', dpi=300)

plt.show()

bundle = None