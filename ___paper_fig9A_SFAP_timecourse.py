import PyPNS
import matplotlib.pyplot as plt
import numpy as np
import os

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 25000 # um

# bundle guide
bundleGuide = PyPNS.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon=30)

# axon definitions
myelinatedParameters = {'fiberD': 3}
unmyelinatedParameters = {'fiberD': 3}

# ----------------------------- stimulation params ---------------------------

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 50., #50,  # Pulse amplitude (mA)
                           'frequency': 20.,  # Frequency of the pulse (kHz)
                           'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                           'stimDur': 0.05,  # Stimulus duration (ms)
                           'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                           'delay': 0.,  # ms
                           # 'invert': True,
                           # 'timeRes': timeRes,
                           }


intraParameters = {'stimulusSignal': PyPNS.signalGeneration.rectangular(**rectangularSignalParams)}

# plotting parameters
(f, axarr) = plt.subplots(1, 2, figsize=(8.4,4))
colors = np.array(((0.,0.,0.), (230., 159., 0.), (86., 180., 233.), (0., 158., 115.)))/255

legends = ['Unmyelinated', 'Myelinated']
recMechLegends = ['homogeneous', 'FEM', 'idealizedCuff']
recMechMarkers = ['o', 'v']

# iterate over axon types
for i in [0,1]:

    # set all properties of the bundle
    bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                        'length': lengthOfBundle,#bundleLengths[i],  # um Axon length
                        'randomDirectionComponent': 0.,
                        'bundleGuide': bundleGuide,

                        'numberOfAxons': 1,  # Number of axons in the bundle
                        'pMyel': i,  # Percentage of myelinated fiber type A
                        'pUnmyel': 1 - i,  # Percentage of unmyelinated fiber type C
                        'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                        'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C
                        'axonCoords': [0, 0],

                        'tStop': 50, # ms
                        'timeRes': 0.0025, #'variable', #

                        'saveV': False,
                        }

    # create the bundle with all properties of axons and recording setup
    bundle = PyPNS.Bundle(**bundleParameters)

    # create the extracellular space models
    LFPMech = []
    LFPMech.append(PyPNS.Extracellular.homogeneous(sigma=1))
    LFPMech.append(PyPNS.Extracellular.precomputedFEM(bundle.bundleCoords))
    LFPMech.append(PyPNS.Extracellular.analytic(bundle.bundleCoords))

    # spiking through a single electrical stimulation
    bundle.add_excitation_mechanism(PyPNS.StimIntra(**intraParameters))

    if i == 1:

        # for myelinated axons place the electrode above the node
        relPosition = 0
        recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                  'radius': 235,
                                  'positionAlongBundle': np.floor(15000. / bundle.axons[0].lengthOneCycle) *
                                                         bundle.axons[0].lengthOneCycle + bundle.axons[0].lengthOneCycle*relPosition,
                                  'numberOfPoles': 1,
                                  'poleDistance': 1000,
                                  'numberOfPoints': 20,
                                  }
        electrodePos = PyPNS.createGeometry.circular_electrode(**recordingParametersNew)

    else:

        recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                  'radius': 235,
                                  'positionAlongBundle': 12000,
                                  'numberOfPoles': 1,
                                  'poleDistance': 1000,
                                  'numberOfPoints': 20,
                                  }
        electrodePos = PyPNS.createGeometry.circular_electrode(**recordingParametersNew)

    # combine electrode position and extracellular field model into a recording mechanism and add to bundle
    modularRecMechs = [0 for ii in range(3)]
    for recMechInd in range(3):
        modularRecMechs[recMechInd] = PyPNS.RecordingMechanism(electrodePos, LFPMech[recMechInd])
        bundle.add_recording_mechanism(modularRecMechs[recMechInd])

    # run the simulation
    bundle.simulate()

    for recMechInd in range(3):

        # load the SFAP of this recording mechanism
        t, SFAPs = bundle.get_SFAPs_from_file(recMechInd)

        # convert from mV to uV
        SFAPs = SFAPs * 1000

        # plot
        axarr[i].plot(t,SFAPs, label=recMechLegends[recMechInd], color=colors[recMechInd+1,:], linewidth=2)
        axarr[i].set_xlabel('time [ms]')

    if i == 0:
        axarr[i].set_ylabel('$V_{ext}$ [$\mu$V]')
        axarr[i].set_title('Unmyelinated')
    else:
        axarr[i].set_title('Myelinated')

plt.legend(loc='best', frameon=False)

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

# zoom in and save figure
# first dimension: zoom, second dimension: un-/myel
xlims = [[(0.8, 38.0), (0.1, 2.5)], [(15.6, 20.7), (0.5, 1.5)]]
ylims = np.array([[(-0.007, 0.0032), (-0.007, 0.0032)], [(-0.0006, 0.00011), (-0.0003, 0.00025)]])*1000

if not os.path.exists('figures'):
    os.makedirs('figures')

figureNames = ['SFAPglobal.eps', 'SFAPzoomed.eps']
for zoomInd in range(2):

    xlimits = xlims[zoomInd]
    ylimits = ylims[zoomInd]

    for axInd in range(2):
        axarr[axInd].set_xlim(xlimits[axInd])
        axarr[axInd].set_ylim(ylimits[axInd])

    plt.tight_layout()
    plt.savefig(os.path.join('figures', figureNames[zoomInd]), format='eps', dpi=300)

plt.show()