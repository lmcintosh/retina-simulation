from numpy import *
from scipy import *
from matplotlib.pyplot import *
import pdb
#from brian import *

def RGCresponse(stimulus, RFradius, RFloc, RGCtype, sampleRate):
    # INPUTS: stimulus x*y*time/sampleRate
    # INPUTS: RFradius in [0,1] space
    # INPUTS: RFloc, RF center in [0,1]x[0,1] space
    # INPUTS: RGCtype, 0 if ON, 1 if OFF
    
    # parameters
    T = stimulus.shape[2]*sampleRate   # duration of stimulus in s
    filterDuration = int(32/(1000/(1/sampleRate))) # must be integer
    nlPoint = -0.9
    nlSlope = 2.
    maxFiringRate = 100.
    #sampleRate = T/stimulus.shape[2]
    onWeights  = .7  # center-surround
    offWeights = -.2    # center-surround
    frequency  = 1.   # of biphasic filter
    phase      = -0.5    # of biphasic filter
    variance   = 1.    # variance of biphasic filter

    # generate stimulus

    #pdb.set_trace()
    # assuming RFradius and RFloc are in [0,1] x [0,1] space
    centerRadius = 0.5*RFradius*mean([stimulus.shape[0],stimulus.shape[1]])
    yLength = stimulus.shape[0]
    xLength = stimulus.shape[1]
    pixelCenter = [int(RFloc[0]*xLength), int(RFloc[1]*yLength)]
    pixelWidth  = [max(int(pixelCenter[0]-centerRadius),int(0)), min(int(pixelCenter[0]+centerRadius),int(stimulus.shape[0]-1))]
    pixelHeight = [max(int(pixelCenter[1]-centerRadius),int(0)), min(int(pixelCenter[1]+centerRadius),int(stimulus.shape[1]-1))]

    if RGCtype == 0:    # ON type
        spatialOutput = zeros(stimulus.shape[2])
        for i in arange(pixelWidth[0], pixelWidth[1]+1):
            for j in arange(pixelHeight[0], pixelHeight[1]+1):
                if ((i-pixelCenter[0])**2+(j-pixelCenter[1])**2) < (centerRadius/2)**2:
                    spatialOutput += stimulus[i,j,:]*onWeights
                elif ((i-pixelCenter[0])**2+(j-pixelCenter[1])**2) < (centerRadius)**2:
                    spatialOutput += stimulus[i,j,:]*offWeights
    elif RGCtype == 1:  # OFF type
        spatialOutput = zeros(stimulus.shape[2])
        for i in arange(pixelWidth[0], pixelWidth[1]+1):
            for j in arange(pixelHeight[0], pixelHeight[1]+1):
                if ((i-pixelCenter[0])**2+(j-pixelCenter[1])**2) < (centerRadius/2)**2:
                    spatialOutput += stimulus[i,j,:]*offWeights
                elif ((i-pixelCenter[0])**2+(j-pixelCenter[1])**2) < (centerRadius)**2:
                    spatialOutput += stimulus[i,j,:]*onWeights
   #       spatialFilter = zeros((stimulus.shape[0],stimulus.shape[1]))
 #       for i in     arange(pixelWidth[0],  pixelWidth[1]-1):
 #           for j in arange(pixelHeight[0], pixelHeight[1]-1):
 #               if sqrt((i-pixelCenter[0])**2+(j-pixelCenter[1])**2) < centerRadius:
 #                   spatialFilter[i,j]     = onWeights
 #                   if sqrt((i-pixelCenter[0])**2+(j-pixelCenter[1])**2) < centerRadius/2:
 #                       spatialFilter[i,j] = offWeights

#    spatialOutput = zeros(stimulus.shape[2])
#    for t in range(stimulus.shape[2]):
#        projection = array(stimulus[:,:,t])*array(spatialFilter)
#        spatialOutput[t] = sum(projection)

    #pdb.set_trace()

    kernel = linearKernel(frequency,phase,variance,filterDuration)
    linearOutput = convolve(spatialOutput,kernel,'valid')   # note this is max - min + 1
    #nonlinearOutput = threshold(linearOutput, nlPoint, nlSlope)

    # normalize nonlinear output to be between 0 and 1
    #nonlinearOutput = nonlinearOutput - min(nonlinearOutput)
    #nonlinearOutput = nonlinearOutput/max(nonlinearOutput)
    nonlinearOutput = exp(linearOutput)/(1+exp(linearOutput))

    spikes = []
    spikeTimes = []
    #timesGen = drange(0,T,sampleRate)
    times = arange(0,len(nonlinearOutput))  # was stimulus.shape[2]
    times = times*sampleRate
    #times = ["%g" % x for x in timesGen]
    print "{} time samples".format(len(times))
    for t in range(len(nonlinearOutput)):
        spikes.append(random.poisson(maxFiringRate*sampleRate*nonlinearOutput[t]))
        if spikes[-1]:
            now = 1000.*array(times[t], float64) + random.rand(1)
            spikeTimes.append(now)

    return spikeTimes

def linearKernel(freq,phase,variance,width):
    #needs to be fixed to be a proper high pass filter
    x = linspace(0,2*pi,width)
    sinusoid = sin(freq*x + phase)
    y = gaussian(0,variance,10e3)
    gauss = histogram(y,width)[0]*1.0/10e3
    kernel = sinusoid*gauss
    kernel = kernel/max(kernel)
    return kernel


def threshold(linearResponse, nlPoint, nlSlope):
    result = zeros(size(linearResponse))

    for i in range(size(linearResponse)):
        if linearResponse[i] >= nlPoint:
            result[i] = nlSlope*(linearResponse[i] - nlPoint)
    return result


def gaussian(mean,var,numSamples):
    samples = sqrt(var)*random.randn(numSamples) + mean
    return samples

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

