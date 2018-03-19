import PyPNS
import PyPNS.analyticFnGen
import matplotlib.pyplot as plt
import numpy as np
import os
import cPickle as pickle

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 5., #50,  # Pulse amplitude (mA)
                           'frequency': 20.,  # Frequency of the pulse (kHz)
                           'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                           'stimDur': 0.05,  # Stimulus duration (ms)
                           'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                           'delay': 0.,  # ms
                           # 'invert': True,
                           # 'timeRes': timeRes,
                           }

intraParameters = {'stimulusSignal': PyPNS.signalGeneration.rectangular(**rectangularSignalParams)}



# axon properites
diametersUnmyel = np.arange(0.2, 1.5, 0.1)
diametersMyel = np.arange(1., 3.5, 0.4)
diametersBothTypes = [diametersUnmyel, diametersMyel]

# plotting
tStartPlots = [0.2, 0.05]

# C-fibres are slower, therefore take more time to reach the electorde
simTimes = [160, 40]
dt = 0.0025

cuffWidthsBothFibres = [np.logspace(-4,-2,19)*2, np.logspace(-4,-1,19)]

# dictionary to save results in
saveDict = {'cuffWidthsUnmyelinated': cuffWidthsBothFibres[0],
            'cuffWidthsMyelinated': cuffWidthsBothFibres[1],
            'unmyelinatedSFAPs': [[] for j in range(len(diametersUnmyel))],
            'myelinatedSFAPs': [[] for j in range(len(diametersMyel))],
            'unmyelinatedDiameters' : diametersUnmyel,
            'myelinatedDiameters': diametersMyel,
            }


legends = ['Unmyelinated', 'Myelinated']
bundleLengths = [25000, 110000]
for i in [0,1]:

    vAPCollection = []

    diameters = diametersBothTypes[i]
    bundleLength = bundleLengths[i]
    cuffWidths = cuffWidthsBothFibres[i]

    # recMechLegends = ['homogeneous', 'FEM', 'idealizedCuff']
    # recMechMarkers = ['o', 'v']

    SFAPsHomo = []
    SFAPsFEM = []
    SFAPsIdeal = []
    CVs = []

    for diameterInd, diameter in enumerate(diameters):

        myelinatedDiam = diameter  # {'distName' : 'normal', 'params' : (1.0, 0.7)}
        unmyelinatedDiam = diameter  # {'distName' : 'normal', 'params' : (0.7, 0.3)}

        # axon definitions
        myelinatedParameters = {'fiberD': myelinatedDiam}
        unmyelinatedParameters = {'fiberD': unmyelinatedDiam}

        # set all properties of the bundle
        bundleParameters = {'radius': 190,  #um radius of the bundle
                            'length': bundleLength,  # um Axon length

                            'numberOfAxons': 1,  # Number of axons in the bundle
                            'pMyel': i,  # Percentage of myelinated fiber type A
                            'pUnmyel': 1 - i,  # Percentage of unmyelinated fiber type C
                            'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                            'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                            'tStop': simTimes[i],
                            'timeRes': dt,

                            'saveV': False,
                            # 'saveLocation': '/path/to/save/directory'
                            }

        # create the bundle with all properties of axons and recording setup
        bundle = PyPNS.Bundle(**bundleParameters)

        # create the extracellular media
        LFPMechs = []
        for cuffWidthInd, cuffWidth in enumerate(cuffWidths):
            LFPMechs.append(PyPNS.Extracellular.analytic(bundle.bundleCoords, interpolator=PyPNS.analyticFnGen.idealizedCuff(cuffWidth)))

        # spiking through a single electrical stimulation
        bundle.add_excitation_mechanism(PyPNS.StimIntra(**intraParameters))

        # recording position
        elecDist = 0
        if i == 1:

            elecDist = np.floor(bundleLength*0.6 / bundle.axons[0].lengthOneCycle) * bundle.axons[0].lengthOneCycle

            recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                      'radius': 235,
                                      'positionAlongBundle': elecDist,
                                      'numberOfPoles': 1,
                                      'poleDistance': 1000,
                                      }

        else:

            elecDist = bundleLength*0.5

            recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                      'radius': 235,
                                      'positionAlongBundle': elecDist,
                                      'numberOfPoles': 1,
                                      'poleDistance': 1000,
                                      }

        electrodePos = PyPNS.createGeometry.circular_electrode(**recordingParametersNew)

        # compose extracellular medium model and recording electrode to a recording mechanism
        modularRecMechs = []
        for recMechIndex in range(len(cuffWidths)):
            modularRecMechs.append(PyPNS.RecordingMechanism(electrodePos, LFPMechs[recMechIndex]))
            bundle.add_recording_mechanism(modularRecMechs[-1])

        # run the simulation
        bundle.simulate()

        SFAPs = []
        for recMechIndex in range(len(cuffWidths)):
            t, SFAP = bundle.get_SFAPs_from_file(recMechIndex)

            SFAPs.append(np.squeeze(SFAP))

        if i == 0:
            saveDict['unmyelinatedSFAPs'][diameterInd] = SFAPs
        else:
            saveDict['myelinatedSFAPs'][diameterInd] = SFAPs

        bundle.clear_all_CAP_files()
        bundle = None

if not os.path.exists('SFAPs'):
    os.makedirs('SFAPs')

pickle.dump(saveDict, open(os.path.join('SFAPs', 'SFAPsCuffWidth3.dict'), "wb"))