import cPickle as pickle
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

colors = np.array(((0.,0.,0.), (230., 159., 0.), (86., 180., 233.), (0., 158., 115.)))/255

# experimentally recorded data for comparison
data = np.loadtxt('experimentalData/pw1ms_amp20.0mA.txt', delimiter=',')

# convert time to ms, voltage mean and standard deviation to mV
timeExp = data[:,0]*1000
voltage = data[:,1]/1000
voltageError = data[:,2]/1000

# load simulated CAP
simOut = np.load('plottedCAP')
tCAP = simOut[0]
CAPs = simOut[1:4]

# time windows for myelinated and unmyelianted signal segments
timesMyel = [8,25]
timesUnmyel = [40, 110]

fieldStrings = ['homogeneous', 'radially inhomgeneous', 'cuff']

# plotting
f, axarr = plt.subplots(2, 2, sharey='row')

# less detail for non-zoomed plots
f_rec_exp_myel, Pxx_den_rec_exp_myel = signal.welch(
    data[np.logical_and(timeExp > timesMyel[0], timeExp < timesMyel[1]), 1] - np.mean(data[:, 1]),
    1. / (timeExp[1] - timeExp[0])*1000, nperseg=32)
f_rec_exp_unmyel, Pxx_den_rec_exp_unmyel = signal.welch(
    data[np.logical_and(timeExp > timesUnmyel[0], timeExp < timesUnmyel[1]), 1] - np.mean(data[:, 1]),
    1. / (timeExp[1] - timeExp[0])*1000, nperseg=32)

# plot experimental recording spectrum
for rowInd in range(2):
    axarr[rowInd][1].plot(f_rec_exp_myel/1000, 20 * np.log10(
        np.sqrt(Pxx_den_rec_exp_myel / np.max(Pxx_den_rec_exp_myel))), label='experiment', color='k', linewidth=2)
    axarr[rowInd][0].plot(f_rec_exp_unmyel/1000, 20 * np.log10(
        np.sqrt(Pxx_den_rec_exp_unmyel / np.max(Pxx_den_rec_exp_unmyel))), label='experiment', color='k', linewidth=2)

# plot simulated spectrum
timeIndMyel = np.logical_and(tCAP > timesMyel[0], tCAP < timesMyel[1])

for fieldTypeInd in range(3):
    CAP = CAPs[fieldTypeInd]
    f_sim_myel, Pxx_sim_myel = signal.welch(CAP[timeIndMyel] - np.mean(CAP[timeIndMyel]), 1000. / (tCAP[1] - tCAP[0]), nperseg=1000)

    timeIndUnmyel = np.logical_and(tCAP > timesUnmyel[0], tCAP < timesUnmyel[1])
    f_sim_unmyel, Pxx_sim_unmyel = signal.welch(CAP[timeIndUnmyel] - np.mean(CAP[timeIndUnmyel]), 1000. / (tCAP[1] - tCAP[0]), nperseg=2000)

    for rowInd in range(2):
        axarr[rowInd][1].plot(f_sim_myel/1000, 20*np.log10(np.sqrt(Pxx_sim_myel / np.max(Pxx_sim_myel))), label=fieldStrings[fieldTypeInd],
                          color=colors[fieldTypeInd + 1], linewidth=2)
        axarr[rowInd][0].plot(f_sim_unmyel/1000, 20*np.log10(np.sqrt(Pxx_sim_unmyel / np.max(Pxx_sim_unmyel))), label=fieldStrings[fieldTypeInd],
                          color=colors[fieldTypeInd + 1], linewidth=2)

axarr[0][0].set_ylim((-55, 5))
axarr[0][0].set_xlim((0, 10))
axarr[0][1].set_xlim((0, 10))
axarr[1][0].set_ylim((-35, 5))
axarr[1][0].set_xlim((0, 2))
axarr[1][1].set_xlim((0, 3.5))

axarr[1][1].legend(loc='best', frameon=False)

axarr[0][0].set_ylabel('normalized power (dB)')
axarr[1][0].set_ylabel('normalized power (dB)')
axarr[1][0].set_xlabel('frequency (kHz)')
axarr[1][1].set_xlabel('frequency (kHz)')

# remove unneccessary axis lines
for ax in axarr.flatten():
    for loc, spine in ax.spines.items():
        if loc in ['right', 'top']:
            spine.set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

if not os.path.exists('figures'):
    os.makedirs('figures')

plt.savefig(os.path.join('figures', 'fig11_spectrum.eps'),
        format='eps', dpi=300)

plt.show()