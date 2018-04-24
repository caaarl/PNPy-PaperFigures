import PyPNS
import matplotlib.pyplot as plt
import numpy as np
import os
import cPickle as pickle

# two simulations, one for myelinated fibres, one for unmyelinated ones.

electrodeDistance = 70000 # um

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------- Simulation parameters are the same for both --------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------- Axon properties ------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# unmyelinated axon diameter distribution
diametersUnmyel = (0.120, 0.17,  0.21,  0.26,  0.32,  0.37,  0.41,  0.47,  0.51,  0.56,  0.62,  0.67,  0.72,  0.77,  0.84,  0.92,  0.97,  1.02,  1.07,  1.12,  1.17,  1.22,  1.27,  1.32, 1.36, 1.41, 1.48, 1.52)
fibreProbabilityUnmyel = (0.0691040631732923, 0.182192465406599, 0.429980837522710, 0.632957475186409, 2.05015339910575,
                              3.10696898591111,  4.54590886074274,  7.22064649366380,  7.60343269800399,  8.61543655035694,
                              8.07683524571988,  7.15617584468796,  7.04457675416097,  6.77590492234067,  5.67583310442061,
                              5.20464797635635,  3.22856301277829,  2.51011904564906,  2.06140597644239,  1.50026642131635,
                              1.32118496258518,  0.849999834520921, 0.760773515404445, 0.312027350382088, 0.200593738933586,
                              0.201222559431810, 0.15, 0.1)

# axon definitions
myelinatedParameters = {'fiberD': {'distName': 'normal', 'params': (1.7, 0.4)}}
unmyelinatedParameters = {'fiberD': {'distName': 'manual', 'params':
    {'diameters': diametersUnmyel, 'densities':fibreProbabilityUnmyel}}}

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------- Load and plot data for reference--------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

data = np.loadtxt('experimentalData/pw1ms_amp20.0mA.txt', delimiter=',')

# convert time to ms, voltage mean and standard deviation to mV
time = data[:,0]*1000
voltage = data[:,1]/1000
voltageError = data[:,2]/1000

# plot the reference
plt.plot(time, voltage, color='black', linewidth=1, label='experimental data')
plt.fill_between(time, voltage-voltageError, voltage+voltageError, facecolor=np.ones(3)*0.6)

# ----------------------------------------------------------------------------------------------------------------------
# --------------------------------------- Unmyelinated simulation-------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# simulation parameters
dt = 0.0025
tStop = 120

# bundle parameters
elecDist = electrodeDistance/5 # unmyelinated fibers conduct too slow, therefore artificially shorten the bundle
nAxons = 2000
bundleLength = elecDist + 12000 # add some length so there is no artefact of the signal reaching the axon end

# set all properties of the bundle
bundleParameters = {'radius': 190, #um radius of the bundle
                    'length': bundleLength,  # um Axon length
                    'randomDirectionComponent': 0.2,

                    'numberOfAxons': nAxons,  # Number of axons in the bundle
                    'pMyel': 0.0,  # Percentage of myelinated fiber type A
                    'pUnmyel': 1.0,  # Percentage of unmyelinated fiber type C
                    'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                    'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                    'tStop': tStop,
                    'timeRes': dt,

                    'saveV': False,
                    # 'saveLocation': '/path/to/save/directory'
                    }

# create the bundle with all properties of axons and recording setup
bundle = PyPNS.Bundle(**bundleParameters)

# create the extracellular media
LFPMech = []
LFPMech.append(PyPNS.Extracellular.homogeneous(sigma=1))
LFPMech.append(PyPNS.Extracellular.precomputedFEM(bundle.bundleCoords))
LFPMech.append(PyPNS.Extracellular.analytic(bundle.bundleCoords))

# spiking through a single electrical stimulation
bundle.add_excitation_mechanism(PyPNS.StimIntra(**intraParameters))


recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                          'radius': 235,
                          'positionAlongBundle': elecDist,
                          'numberOfPoles': 2,
                          'poleDistance': 3000,
                          }

electrodePos = PyPNS.createGeometry.circular_electrode(**recordingParametersNew)

# compose extracellular medium model and recording electrode to a recording mechanism
modularRecMechs = []
for recMechIndex in range(3):
    modularRecMechs.append(PyPNS.RecordingMechanism(electrodePos, LFPMech[recMechIndex]))
    bundle.add_recording_mechanism(modularRecMechs[-1])

# run the simulation
bundle.simulate()

CAPUnmyelinated = [[] for i in range(3)]
for recMechIndex in range(3):
    timeUnmyelinated, CAPUnmyelinated[recMechIndex] = bundle.get_CAP_from_file(recMechIndex)

bundle = None

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------- Myelinated simulation-------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# simulation parameters
dt = 0.0025
tStop = 150

# bundle parameters
elecDist = electrodeDistance # myelinated fibers conduct too slow, therefore artificially shorten the bundle
nAxons = 150
bundleLength = elecDist + 12000 # add some length so there is no artefact of the signal reaching the axon end

# set all properties of the bundle
bundleParameters = {'radius': 190, #um radius of the bundle
                    'length': bundleLength,  # um Axon length
                    'randomDirectionComponent': 0.2,

                    'numberOfAxons': nAxons,  # Number of axons in the bundle
                    'pMyel': 1.0,  # Percentage of myelinated fiber type A
                    'pUnmyel': 0.0,  # Percentage of unmyelinated fiber type C
                    'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                    'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                    'tStop': tStop,
                    'timeRes': dt,

                    'saveV': False,
                    # 'saveLocation': '/media/.../PyPNS'
                    }

# create the bundle with all properties of axons and recording setup
bundle = PyPNS.Bundle(**bundleParameters)

# create the extracellular media
LFPMech = []
LFPMech.append(PyPNS.Extracellular.homogeneous(sigma=1))
LFPMech.append(PyPNS.Extracellular.precomputedFEM(bundle.bundleCoords))
LFPMech.append(PyPNS.Extracellular.analytic(bundle.bundleCoords))

# spiking through a single electrical stimulation
bundle.add_excitation_mechanism(PyPNS.StimIntra(**intraParameters))


recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                          'radius': 235,
                          'positionAlongBundle': elecDist,
                          'numberOfPoles': 2,
                          'poleDistance': 3000,
                          }

electrodePos = PyPNS.createGeometry.circular_electrode(**recordingParametersNew)

# compose extracellular medium model and recording electrode to a recording mechanism
modularRecMechs = []
for recMechIndex in range(3):
    modularRecMechs.append(PyPNS.RecordingMechanism(electrodePos, LFPMech[recMechIndex]))
    bundle.add_recording_mechanism(modularRecMechs[-1])

# run the simulation
bundle.simulate()

# save the CAP
CAPMyelinated = [[] for i in range(3)]
for recMechIndex in range(3):
    timeMyelinated, CAPMyelinated[recMechIndex] = bundle.get_CAP_from_file(recMechIndex)

bundle = None

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------- Combine CAPs and plot ------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

legends = ['Unmyelinated', 'Myelinated']
recMechLegends = ['homogeneous', 'radially inhomogeneous', 'cuff']
recMechMarkers = ['o', 'v']

# number of time steps per fibre type
nUnmyelinated = np.shape(timeUnmyelinated)[0]
nMyelinated = np.shape(timeMyelinated)[0]
# max number of time steps
nMax = np.max((nMyelinated, nUnmyelinated))

if nUnmyelinated >= nMyelinated:
    t = timeUnmyelinated
else:
    t = timeMyelinated


CAPs = []
for recMechIndex in range(3):

    # array to accommodate both waveforms
    CAP = np.zeros(nMax)

    CAP[0:nMyelinated] = CAPMyelinated[recMechIndex]

    CAP[0:nUnmyelinated] += CAPUnmyelinated[recMechIndex]

    plt.plot(t, CAP, label=recMechLegends[recMechIndex])

    CAPs.append(CAP)


plt.legend(loc='best', frameon=False)
plt.xlabel('time (ms)')
plt.ylabel('extracellular voltage (mV)')

if not os.path.exists('figures'):
    os.makedirs('figures')

plt.savefig(os.path.join('figures', 'fig10_CAP.eps'),
        format='eps', dpi=300)

np.save(open('plottedCAP', 'wb'), np.row_stack((t, CAPs)))

plt.show()
