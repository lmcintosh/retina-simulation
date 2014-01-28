import numpy
#import brian
import retinaFunctions
import pdb


class rf:
    def __init__(self, center, radius):
        self.center = center
        self.radius = radius

    def overlap_area(self, other_rf):
        dist = numpy.linalg.norm(self.center-other_rf.center)
        if dist > self.radius+other_rf.radius:
            return 0.0
        elif dist + self.radius <= other_rf.radius:
            return 1.0*numpy.pi*self.radius**2
        elif dist + other_rf.radius <= self.radius:
            return 1.0*numpy.pi*other_rf.radius**2
        else:
            jointtheta = numpy.arccos((self.radius**2+other_rf.radius**2-dist**2)/(2*self.radius*other_rf.radius))
            selftheta = numpy.arcsin(self.radius*numpy.sin(jointtheta)/dist) 
            othertheta = numpy.arcsin(other_rf.radius*numpy.sin(jointtheta)/dist) 
            
            h = self.radius * numpy.sin(selftheta)
            total_area = dist*h/2
            overlapped_area = (selftheta/(2*numpy.pi))*numpy.pi*self.radius**2 + (othertheta/(2*numpy.pi))*numpy.pi*other_rf.radius**2 
            #selfadj = self.radius * numpy.cos(selftheta)
            #otheradj = dist - selfadj

            return 2.0*(overlapped_area - total_area)

    def find_overlapping_indices(self, other_rf_list, prev_min_row = 0, prev_min_col = 0):
        indices = []
        lookback = 4
        row = max(prev_min_row - lookback, 0)
        while row < len(other_rf_list):
            col = max(prev_min_col - lookback, 0)
            while col < len(other_rf_list[row]):
                curr_overlap = self.overlap_area(other_rf_list[row][col])
                if curr_overlap > 0:
                    indices.append([row,col,curr_overlap])
                elif len(indices)>0 and indices[-1][0] == row and indices[-1][1] < col - lookback:
                    break
                col += 1
            if len(indices)>0 and indices[-1][0] < row - lookback:
                break
            row += 1

        return indices


class connection:
    def __init__(self, syn_type, syn_weight, from_id, to_id):
        self.syn_type = syn_type
        self.syn_weight = syn_weight
        self.from_id = from_id
        self.to_id = to_id

class v1_layer:
    def __init__(self,rf_type='simple3',orientation=0):
        self.inputs = []
        self.outputs = []

class circular_rf_layer:
    def __init__(self, rf_xlim,rf_ylim,width,height,std_offset=[0.025], mean_radius = [0.1], std_radius = [0.0001]):
        self.rf_xlim = rf_xlim
        self.rf_ylim = rf_ylim
        self.width = width
        self.height = height
        self.std_offset = std_offset
        self.mean_radius = mean_radius
        self.std_radius = std_radius
        self.rfs = [[[]]]
        self.ids = [[[]]]

    def tessellate(self,id_start):
        cell_id = id_start
        self.rfs = []
        self.ids = []
        for type_index in range(len(self.std_offset)):
            x_locs = numpy.linspace(self.rf_xlim[0],self.rf_xlim[1],self.width)
            y_locs = numpy.linspace(self.rf_ylim[0],self.rf_ylim[1],self.height)
            #x_locs = numpy.arange(self.rf_xlim[0],self.rf_xlim[1],(1-self.overlap[type_index])*self.mean_radius[type_index]*2)
            #y_locs = numpy.arange(self.rf_ylim[0],self.rf_ylim[1],(1-self.overlap[type_index])*self.mean_radius[type_index]*numpy.sqrt(3))
            k = self.mean_radius[type_index] ** 2 / self.std_radius[type_index]
            theta = self.std_radius[type_index] / self.mean_radius[type_index]
            self.rfs.append([])
            self.ids.append([])
            for i in range(len(y_locs)):
                self.rfs[type_index].append([])
                self.ids[type_index].append([])
                for j in range(len(x_locs)):
                    location = numpy.array((x_locs[j]+(i%2)*self.mean_radius[type_index],y_locs[i]))
                    rand_location = location + numpy.random.randn(2) *self.std_offset[type_index]
                    radius = numpy.random.gamma(k,theta)
                    self.rfs[type_index][i].append(rf(rand_location,radius))
                    self.ids[type_index][i].append(cell_id)
                    cell_id += 1
        self.id_range = [id_start,cell_id]
        return cell_id

class RGC_layer:
    def __init__(self, rf_xlim = numpy.array([0,1]),rf_ylim = numpy.array([0,1]),types=[0,1,2,3],height=8,width=8):
        self.rf_xlim = rf_xlim
        self.width = width
        self.rf_ylim = rf_ylim
        self.height = height
        self.types = types
        self.rf_layers  = circular_rf_layer(self.rf_xlim,self.rf_ylim,self.width,self.height,std_offset=[0.0025, 0.0025, 0.0075, 0.0075],mean_radius=[1./(width),1./(width),3./(width),3./(width)],std_radius=[0.00001,0.00001,0.00003,0.00003])

    def construct_layer(self,id_start=0):
        return self.rf_layers.tessellate(id_start)

    def generate_spikes(self,stimulus,sampleRate):
        spikes = []
        for t in range(len(self.types)):
            for row in range(len(self.rf_layers.rfs[t])):
                for col in range(len(self.rf_layers.rfs[t][row])):
                    current_RGC = self.rf_layers.rfs[t][row][col]
                    #pdb.set_trace()
                    spikes.append([t, row, col, retinaFunctions.RGCresponse(stimulus,current_RGC.radius,current_RGC.center,t%2,sampleRate)])
        return spikes


#class LGN_layer:
#    def __init__(self, rf_xlim = numpy.array([0,1]),rf_ylim = numpy.array([0,1]),lagged_percentage = [0.4, 0.4, 0.05, 0.05],types=[4,5,6,7]):
#        self.rf_xlim = rf_xlim
#        self.rf_ylim = rf_ylim
#        self.types = types
#        self.rf_layers  = circular_rf_layer(self.rf_xlim,self.rf_ylim,std_offset=[0.025, 0.025, 0.075, 0.075],overlap=[0.25,0.25,0.75,0.75],mean_radius=[0.1,0.1,0.3,0.3],std_radius=[0.0001,0.0001,0.0003,0.0003])
#
#    def construct_layer(self,id_start=0):
#        return self.rf_layers.tessellate(id_start)
#
#    def connect_RGC(self,RGC):
#        #for each type of RGC
#        connections = []
#        for type_index in range(len(self.types)):
#            prev_min_row = 0
#            for row in range(len(self.rf_layers.rfs[type_index])):
#                prev_min_col = 0
#                for col in range(len(self.rf_layers.rfs[type_index][row])):
#                    LGN_rf = self.rf_layers.rfs[type_index][row][col]
#                    RGC_rfs = RGC.rf_layers.rfs[type_index]
#                    indices = numpy.array(LGN_rf.find_overlapping_indices(RGC_rfs,prev_min_row,prev_min_col))
#                    if len(indices) == 0:
#                        print 'no overlapping indices found type {} row {} col {}'.format(type_index,row,col)
#                        continue
#
#                    max_element = numpy.argmax(indices[:,2]) 
#                    max_row = int(indices[max_element,0])
#                    max_col = int(indices[max_element,1])
#                    connections.append(connection(0,1,RGC.rf_layers.ids[type_index][max_row][max_col],self.rf_layers.ids[type_index][row][col]))
#
#                    prev_min_row = numpy.argmin(indices[:,0]) 
#                    prev_min_col = numpy.argmin(indices[:,1]) 
#        return connections

#class PGN_layer:
#    def __init__(self, rf_xlim = numpy.array([-1,1]),rf_ylim = numpy.array([-1,1]):
#        self.rf_xlim = rf_xlim
#        self.rf_ylim = rf_ylim
#        self.p_on  = circular_rf_layer(self.rf_xlim,self.rf_ylim,std_offset=0.025,mean_radius=0.1,std_radius=0.0001)
#        self.p_off = circular_rf_layer(self.rf_xlim,self.rf_ylim,std_offset=0.025,mean_radius=0.1,std_radius=0.0001)
#        self.m_on  = circular_rf_layer(self.rf_xlim,self.rf_ylim,std_offset=0.075,mean_radius=0.3,std_radius=0.0003)
#        self.m_off = circular_rf_layer(self.rf_xlim,self.rf_ylim,std_offset=0.075,mean_radius=0.3,std_radius=0.0003)


