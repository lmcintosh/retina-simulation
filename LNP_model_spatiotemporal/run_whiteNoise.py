import visual_neuron
import pylab
import numpy
import pdb
import matplotlib
import retinaFunctions
import PIL
import Image
from scipy import *

numSamples = 5000;
sampleRate = 1./500;

#stimulus = numpy.zeros((100,100,50))
random.seed(42)

side_length = 100
stimulus = numpy.zeros((side_length,side_length,numSamples))
for x in range(0,side_length/2,2):
    stimulus[2*x:2*x+2,:,:] = 1
pylab.imshow(stimulus[:,:,0])
pylab.show()
    
#pic = Image.open("./Natural_Images/DSC_0001.JPG")
#
# crop image
#left   = 0
#top    = 0
#width  = 1000
#height = 1000
#box    = (left, top, left+width, top+height)
#pic    = pic.crop(box)

#pixRGB = numpy.array(pic)
#pixGS = 0.299*pixRGB[:,:,0] + 0.587*pixRGB[:,:,1] + 0.114*pixRGB[:,:,2]
#for i in range(pixRGB.shape[0]):
#    for j in range(pixRGB.shape[1]):
#        pixGS[i,j] = int(0.299*pixRGB[i,j,0] + 0.587*pixRGB[i,j,1] + 0.114*pixRGB[i,j,2])
#print pixGS.shape

#stimulus = zeros([pixGS.shape[0], pixGS.shape[1], numSamples])
#stimulus = numpy.tile(pixGS,[1,1,numSamples])
#for z in range(stimulus.shape[2]):
#    stimulus[:,:,z] = pixGS
height = 16 
width = 16 

RGC = visual_neuron.RGC_layer(width=width,height=height) # specify RF, type
print 'constructing RGC'
stop_id = RGC.construct_layer(id_start = 0)
print 'generating RGC spikes'
responses = RGC.generate_spikes(stimulus,sampleRate)

# save each RGC type, RF center (int space), and spiketrain time stamps
# into a new file
chip_height = 65
chip_width = 64
row_offset = int(chip_height/4-height/2)
col_offset = int(chip_width/4-width/2)
for i in range(len(responses)):
    tp = responses[i][0]
    row = responses[i][1]
    col = responses[i][2]
    spikes = responses[i][3]
    dataFile = open(str(int(tp)) + '-' + str(int(row+row_offset)) + '-' + str(int(col+col_offset)) + '-' + 'rgc' + '.dat','a')
    for j in range(len(spikes)):
        dataFile.write(str(spikes[j]) + '\n')
        #print spikes[j]
    dataFile.close()
    #print responses[i]



#LGN = visual_neuron.LGN_layer()
#print 'constructing LGN'
#LGN.construct_layer(id_start = stop_id)
#print 'connecting RGC to LGN'
#connections = LGN.connect_RGC(RGC)


#spikes = retinaFunctions.RGCresponse(stimulus, 0.1, [0,0], 0)


#Pe = RGC.generate_brian()

#
#centers_P_ON = numpy.array(test_P.centers)
#centers_P_ON = centers_P.reshape(centers_P.size/2,2)
#radii_P_ON = numpy.array(test_P.radii)
#radii_P_ON = radii_P.reshape(radii_P.size)
#centers_P_OFF = numpy.array(test_P.centers)
#centers_P_OFF = centers_P.reshape(centers_P.size/2,2)
#radii_P_OFF = numpy.array(test_P.radii)
#radii_P_OFF = radii_P.reshape(radii_P.size)
#
#centers_M = numpy.array(test_M.centers)
#centers_M = centers_M.reshape(centers_M.size/2,2)
#radii_M = numpy.array(test_M.radii)
#radii_M = radii_M.reshape(radii_M.size)
#pylab.figure()
#ax = pylab.subplot(221,aspect=1)
#for i in range(centers_P.shape[0]):
#    c = matplotlib.patches.Circle(centers_P[i],radii_P[i],lw=1,ec='b',fill=False)
#    ax.add_patch(c)
#pylab.xlim(test_P.rf_xlim)
#pylab.ylim(test_P.rf_ylim)
#ax = pylab.subplot(222,aspect=1)
#for i in range(centers_M.shape[0]):
#    c = matplotlib.patches.Circle(centers_M[i],radii_M[i],lw=1,ec='r',fill=False)
#    ax.add_patch(c)
#pylab.xlim(test_M.rf_xlim)
#pylab.ylim(test_M.rf_ylim)
#pylab.show()
#
#pdb.set_trace()
