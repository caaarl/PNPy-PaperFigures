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

simTimes = [30,25]

saveDict = {'unmyelinatedDiameters' : diametersUnmyel,
            'myelinatedDiameters' : diametersMyel,
            'unmyelinatedIMem' : [],
            'myelinatedIMem' : []
            }

amplitudes = np.zeros((2,3,len(diametersUnmyel)))

legends = ['Unmyelinated', 'Myelinated']
bundleLengths = [5000, 25000]
for i in [0,1]:

    diameters = diametersBothTypes[i]

    ImemInts = []

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
                            'saveI': True
                            }

        # create the bundle with all properties of axons and recording setup
        bundle = PyPNS.Bundle(**bundleParameters)

        # spiking through a single electrical stimulation
        bundle.add_excitation_mechanism(PyPNS.StimIntra(**intraParameters))

        # run the simulation
        bundle.simulate()

        # load current
        t, imem = bundle.get_imem_from_file_axonwise(0)

        imem = np.array(imem)

        if i == 1:

            nodeIndices = np.where(np.array(bundle.axons[0].segmentTypes) == 'n')[0]
            nodeIndex = nodeIndices[int(len(nodeIndices) * 0.5)]

            current = np.sum(imem[nodeIndex:nodeIndex+3, :], axis=0)
            segmentLength = bundle.axons[0].lengthOneCycle

        else:
            nodeIndex = 10

            current = imem[nodeIndex, :]
            segmentLength = bundle.axons[0].length[0]

        iMemInt = np.trapz(np.abs(current), t) / segmentLength

        ImemInts.append(iMemInt)

        # clear
        bundle = None

    # save everything to a save dictionary
    if i == 0:
        saveDict['unmyelinatedIMem'] = ImemInts
    else:
        saveDict['myelinatedIMem'] = ImemInts

if not os.path.exists('iMems'):
    os.makedirs('iMems')

pickle.dump(saveDict, open(os.path.join('iMems', 'iMemIntVsDiam3.dict'), "wb"))