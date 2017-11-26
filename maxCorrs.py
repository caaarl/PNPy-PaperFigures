import numpy as np


def max_corr (a, v):
    a = (a - np.mean(a)) / (np.std(a) * len(a))
    v = (v - np.mean(v)) / np.std(v)
    return np.max(np.correlate(a, v, "full"))


def max_corrs(signalArray):

    signalArray = np.array(signalArray)

    signalsShape = np.shape(signalArray)
    nSig = min(signalsShape)
    nSamples = max(signalsShape)
    if signalsShape[0] > signalsShape[1]:
        signalArray = np.transpose(signalArray)

    corrs = []

    for i in range(nSig):
        for j in range(i+1, nSig):
            corrs.append(max_corr(signalArray[i,:], signalArray[j,:]))

    return corrs


if __name__ == '__main__':

    nSig = 10
    nSamples = 10000

    # signals = np.random.rand(nSig, nSamples)
    signals = np.ones((nSig, nSamples))

    corrs = max_corrs(signals)

    print 'hm'