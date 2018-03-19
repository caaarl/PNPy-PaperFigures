import PyPNS
import matplotlib.pyplot as plt
import numpy as np
import os

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 7000
numberOfAxons = 1

# bundle guide
segmentLengthAxon = 30
bundleGuide = PyPNS.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

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

# ----------------------------- recording params -------------------------------

recordingParametersNew = {'bundleGuide': bundleGuide,
                          'radius': 250,
                          'positionAlongBundle': 4000,
                          'numberOfPoles': 1,
                          'poleDistance': 1000,
                        }

# random direction components
RDCs = [0.2, 0.6, 1.0]


for RDCInd, RDC in enumerate(RDCs):

    # axon definitions
    myelinatedDiam = 1.
    unmyelinatedDiam = 1.
    myelinatedParameters = {'fiberD': myelinatedDiam}
    unmyelinatedParameters = {'fiberD': unmyelinatedDiam}

    # set all properties of the bundle
    bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                        'length': 7000,  # um Axon length
                        'randomDirectionComponent': RDC,

                        'numberOfAxons': 10,  # Number of axons in the bundle
                        'pMyel': 0.,  # Percentage of myelinated fiber type A
                        'pUnmyel': 1.,  # Percentage of unmyelinated fiber type C
                        'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                        'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                        'tStop': 10,
                        'timeRes': 0.0025, #

                        'saveV': False,
                        }

    # create the bundle with all properties of axons and recording setup
    bundle = PyPNS.Bundle(**bundleParameters)

    ax = PyPNS.plot.geometry_definition(bundle)
    ax.set_xlim((1000, 1300))
    ax.set_ylim((-200, 200))
    ax.set_zlim((-200, 200))
    ax.view_init(90, 90)
    # plt.show()

plt.show()

bundle = None