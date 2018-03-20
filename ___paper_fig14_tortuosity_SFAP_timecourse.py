import PyPNS
import matplotlib.pyplot as plt
import numpy as np
import os
import maxCorrs

# ----------------------------- bundle params -------------------------------

nAxons = 10

# set length of bundle and number of axons
lengthOfBundle = 25000 # um

# bundle guide
bundleGuide = PyPNS.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon=30)

# axon definitions
myelinatedParameters = {'fiberD': 3}
unmyelinatedParameters = {'fiberD': 3}

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
colors = np.array(((0.,0.,0.), (230., 159., 0.), (86., 180., 233.), (0., 158., 115.)))/255

legends = ['Unmyelinated', 'Myelinated']
recMechLegends = ['homogeneous', 'FEM', 'idealizedCuff']
recMechMarkers = ['o', 'v']

RDCs = [0.2, 0.6, 1.0]
tSims = [85, 10]

# iterate over axon types
axs = []
axref = [[] for k in range(3)]

for i in [0,1]:

    maxCorrMean = np.zeros((3, len(RDCs)))
    maxCorrStd = np.zeros((3, len(RDCs)))

    for RDCInd, RDC in enumerate(RDCs):

        # set all properties of the bundle
        bundleParameters = {'radius': 190,  #um Radius of the bundle (typically 0.5-1.5mm)
                            'length': lengthOfBundle, # um Axon length
                            'randomDirectionComponent': RDC,
                            'bundleGuide': bundleGuide,

                            'numberOfAxons': nAxons,  # Number of axons in the bundle
                            'pMyel': i,  # Percentage of myelinated fiber type A
                            'pUnmyel': 1 - i,  # Percentage of unmyelinated fiber type C
                            'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                            'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C
                            'axonCoords': [0, 0],

                            'tStop': tSims[i], # ms
                            'timeRes': 0.0025, # ms

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

        # calculate shape similarity
        for recMechInd in range(3):

            # load the SFAP of this recording mechanism
            t, SFAPs = bundle.get_SFAPs_from_file(recMechInd)

            # compute maximum cross correlation between signals
            corrs = maxCorrs.max_corrs(SFAPs)

            maxCorrMean[recMechInd, RDCInd] = np.mean(corrs)
            maxCorrStd[recMechInd, RDCInd] = np.std(corrs)

        # plotting
        for recMechInd in range(1,3):

            # load the SFAP of this recording mechanism
            t, SFAPs = bundle.get_SFAPs_from_file(recMechInd)

            # convert from mV to uV
            SFAPs = np.array(SFAPs) * 1000

            if RDCInd == 0:
                ax = plt.subplot2grid((4, 6), (recMechInd - 1, RDCInd + i * 3))
                axref[recMechInd] = ax
            else:
                ax = plt.subplot2grid((4,6),(recMechInd-1,RDCInd + i*3), sharey = axref[recMechInd])
            axs.append(ax)

            # plot
            ax.plot(t,SFAPs, label=recMechLegends[recMechInd], color=colors[recMechInd+1,:], linewidth=2)

            if i == 0:
                if recMechInd == 1:
                    if RDCInd == 0:
                        ax.set_xlim((12, 23))
                    elif RDCInd == 1:
                        ax.set_xlim((12, 26))
                    elif RDCInd == 2:
                        ax.set_xlim((12, 80))
                elif recMechInd == 2:
                    if RDCInd == 0:
                        ax.set_xlim((1, 40))
                    elif RDCInd == 1:
                        ax.set_xlim((1, 40))
                    elif RDCInd == 2:
                        ax.set_xlim((1, 80))
            elif i==1:
                if recMechInd == 1:
                    if RDCInd == 0:
                        ax.set_xlim((0.5, 2))
                    elif RDCInd == 1:
                        ax.set_xlim((0.5, 2))
                    elif RDCInd == 2:
                        ax.set_xlim((1, 4.5))
                elif recMechInd == 2:
                    if RDCInd == 0:
                        ax.set_xlim((0, 2.5))
                    elif RDCInd == 1:
                        ax.set_xlim((0, 2.5))
                    elif RDCInd == 2:
                        ax.set_xlim((0, 6))

            if i == 0 and RDCInd == 0:
                ax.set_ylabel('$V_{ext}$ [$\mu$V]')

            if RDCInd == 1:
                ax.set_xlabel('time [ms]')

            if RDCInd > 0:
                plt.setp(ax.get_yticklabels(), visible=False)

    # plot maxcorrs
    ax = plt.subplot2grid((4, 6), (2, 3*i), colspan=3, rowspan=2)
    axs.append(ax)

    if i == 0:
        ax.set_ylabel('pairwise correlation')
    for recMechInd in range(3):
        print 'std'
        print maxCorrStd[recMechInd,:]
        print 'mean'
        print maxCorrMean[recMechInd, :]
        ax.errorbar(RDCs, maxCorrMean[recMechInd,:], yerr=maxCorrStd[recMechInd,:], color=colors[recMechInd+1,:], ecolor=colors[recMechInd+1,:])


# remove unneccessary axis lines
for ax in axs:
    for loc, spine in ax.spines.items():
        if loc in ['right', 'top']:
            spine.set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

plt.tight_layout()

if not os.path.exists('figures'):
    os.makedirs('figures')

plt.savefig(os.path.join('figures', 'fig14.eps'), format='eps', dpi=300)

plt.show()