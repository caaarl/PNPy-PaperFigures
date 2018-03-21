import PyPNS
import matplotlib.pyplot as plt
import numpy as np
import os
import cPickle as pickle

# stimulation parameters
# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 10., #50,  # Pulse amplitude (mA)
                           'frequency': 20.,  # Frequency of the pulse (kHz)
                           'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                           'stimDur': 0.05,  # Stimulus duration (ms)
                           'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                           'delay': 0.,  # ms
                           # 'invert': True,
                           # 'timeRes': timeRes,
                           }
#
intraParameters = {'stimulusSignal': PyPNS.signalGeneration.rectangular(**rectangularSignalParams)}



# ------------------------------------------------------------------------------
# ---------------------------------- CALCULATION -------------------------------
# ------------------------------------------------------------------------------

diametersUnmyel = np.arange(0.2, 4, 0.3)
diametersMyel = diametersUnmyel
diametersBothTypes = [diametersUnmyel, diametersMyel]

simTimes = [100, 20]

saveDict = {'unmyelinatedDiameters' : diametersUnmyel,
            'unmyelinatedSFAPsHomo': [],
            'unmyelinatedSFAPsFEM': [],
            'unmyelinatedCV' : [],
            't': [],
            'myelinatedDiameters': diametersMyel,
            'myelinatedSFAPsHomo': [],
            'myelinatedSFAPsFEM': [],
            'myelinatedCV' : [],
            }

amplitudes = np.zeros((2,3,len(diametersUnmyel)))

legends = ['Unmyelinated', 'Myelinated']
bundleLengths = [25000, 25000]
for i in [0,1]:

    diameters = diametersBothTypes[i]

    SFAPsHomo = []
    SFAPsFEM = []
    SFAPsIdeal = []

    for diameterInd, diameter in enumerate(diameters):

        myelinatedDiam = diameter
        unmyelinatedDiam = diameter

        # axon definitions
        myelinatedParameters = {'fiberD': myelinatedDiam}
        unmyelinatedParameters = {'fiberD': unmyelinatedDiam}

        # set all properties of the bundle
        bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                            'length': bundleLengths[i],  # um Axon length

                            'numberOfAxons': 1,  # Number of axons in the bundle
                            'pMyel': i,  # Percentage of myelinated fiber type A
                            'pUnmyel': 1 - i,  # Percentage of unmyelinated fiber type C
                            'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                            'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                            'tStop': simTimes[i],
                            'timeRes': 0.0025, #'variable', #

                            'saveV': False,
                            }

        # create the bundle with all properties of axons and recording setup
        bundle = PyPNS.Bundle(**bundleParameters)

        # extracellular medium models
        LFPMech = []
        LFPMech.append(PyPNS.Extracellular.homogeneous(sigma=1))
        LFPMech.append(PyPNS.Extracellular.precomputedFEM(bundle.bundleCoords))
        LFPMech.append(PyPNS.Extracellular.analytic(bundle.bundleCoords))

        # spiking through a single electrical stimulation
        bundle.add_excitation_mechanism(PyPNS.StimIntra(**intraParameters))

        # position on node for myelinated fibers, don't care about position for unmyelinated ones
        if i == 1:

            recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                      'radius': 200,
                                      'positionAlongBundle': np.floor(bundleLengths[i]*0.8 / bundle.axons[0].lengthOneCycle) *
                                                             bundle.axons[0].lengthOneCycle,
                                      'numberOfPoles': 1,
                                      'poleDistance': 1000,
                                      }

        else:
            recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                      'radius': 100,
                                      'positionAlongBundle': 12000,
                                      'numberOfPoles': 1,
                                      'poleDistance': 1000,
                                      }

        electrodePos = PyPNS.createGeometry.circular_electrode(**recordingParametersNew)

        # combine electrode position and extracellular medium model
        modularRecMechs = []
        for recMechIndex in [0,1,2]:
            modularRecMechs.append(PyPNS.RecordingMechanism(electrodePos, LFPMech[recMechIndex]))
            bundle.add_recording_mechanism(modularRecMechs[-1])

        # run the simulation
        bundle.simulate()

        # obtain SFAPs and add
        t, SFAPs = bundle.get_SFAPs_from_file(0)
        SFAPNoArt = SFAPs[100:]
        amplitudes[i, 0, diameterInd] = np.max(SFAPNoArt) - np.min(SFAPNoArt)
        SFAPsHomo.append(np.squeeze(SFAPs))

        _, SFAPs = bundle.get_SFAPs_from_file(1)
        SFAPNoArt = SFAPs[100:]
        amplitudes[i, 1, diameterInd] = np.max(SFAPNoArt) - np.min(SFAPNoArt)
        SFAPsFEM.append(np.squeeze(SFAPs))

        _, SFAPs = bundle.get_SFAPs_from_file(2)
        SFAPNoArt = SFAPs[100:]
        amplitudes[i, 2, diameterInd] = np.max(SFAPNoArt) - np.min(SFAPNoArt)
        SFAPsIdeal.append(np.squeeze(SFAPs))

        # clear for memory reasons
        bundle.clear_all_CAP_files()
        bundle = None

    # save everything to a save dictionary
    if i == 0:
        saveDict['unmyelinatedSFAPsHomo'] = SFAPsHomo
        saveDict['unmyelinatedSFAPsFEM'] = SFAPsFEM
        saveDict['unmyelinatedSFAPIdeal'] = SFAPsIdeal
    else:
        saveDict['myelinatedSFAPsHomo'] = SFAPsHomo
        saveDict['myelinatedSFAPsFEM'] = SFAPsFEM
        saveDict['myelinatedSFAPIdeal'] = SFAPsIdeal

if not os.path.exists('SFAPs'):
    os.makedirs('SFAPs')

pickle.dump(saveDict, open(os.path.join('SFAPs', 'SFAPsForAmplitudePlot.dict'), "wb"))