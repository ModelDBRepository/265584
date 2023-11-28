from neuron import h

class Synapse_py3:
    def __init__(self,source,target,section,weight = 1):
		
        self.input = h.NetStim(0.5)
        self.input.start = -10
        self.input.number = 1
        self.input.interval = 1e9
        self.weight = weight


        self.postsyns = {}

        if (type(source) == type('s')):
            sourcetype = source
           
                             
#granule cell     
        if sourcetype == 'mossy':
            if target.whatami == 'GrC_2020_regular':
                self.whatami = "syn_mossytoGrC_det"
                self.postsyns['AMPA'] = [h.GRANULE_Ampa_det_vi(0.9, sec=section)]
                self.postsyns['AMPA'][0].tau_facil=5
                self.postsyns['AMPA'][0].tau_rec=8 #25 
                self.postsyns['AMPA'][0].tau_1=1
                self.postsyns['AMPA'][0].gmax = 1200
                self.postsyns['AMPA'][0].U=0.43

                self.nc_syn = [h.NetCon(self.input,receptor[0],0,0.1,1) for receptor in self.postsyns.values()]
            
        elif sourcetype == 'mossynmda':
            if target.whatami == 'GrC_2020_regular':
                self.whatami = "syn_mossytoGrC_det_nmda"
                self.postsyns['NMDA'] = [h.GRANULE_Nmda_det_vi(0.9, sec=section)]
                self.postsyns['NMDA'][0].tau_facil=5
                self.postsyns['NMDA'][0].tau_rec=8 
                self.postsyns['NMDA'][0].tau_1=1
                self.postsyns['NMDA'][0].gmax = 18800  
                self.postsyns['NMDA'][0].U=0.43
                self.nc_syn = [h.NetCon(self.input,receptor[0],0,0.1,1) for receptor in self.postsyns.values()]
		
        else:
            print('SOURCE TYPE DOES NOT EXIST SOMETHING WRONG!!!!!!!!!')
            
            

    
