import matplotlib as mpl
mpl.use('tkagg')  

from neuron import h,gui 
from Cereb_GrC_regular import Grc_regular
import multiprocessing
import numpy as np
import matplotlib.pyplot as plt

cell = Grc_regular(1)

#h.nrncontrolmenu()

time_step = h.CVode()
time_step.active(0) #0 fixed step, 1 variable time step

cpu = multiprocessing.cpu_count()
h.load_file("parcom.hoc")
p = h.ParallelComputeTool()
if cpu > 8:
    p.change_nthread(8,1)
    print('Maximum 8 threads')
else:
    p.change_nthread(cpu,1)    
    print('N° of treads', cpu)
p.multisplit(1)

h.dt = 0.025
h.celsius = 32
h.tstop = 5000
h.v_init = -70

stim = [h.IClamp(0.5,sec=cell.soma[0]), h.IClamp(0.5,sec=cell.soma[0]), h.IClamp(0.5,sec=cell.soma[0])]

stim[0].delay = 100
stim[0].dur = 1500
stim[0].amp = 0.01 #10pA

stim[1].delay = 1700
stim[1].dur = 1500
stim[1].amp = 0.016 #16pA

stim[2].delay = 3300
stim[2].dur = 1500
stim[2].amp = 0.022 #22pA

h('load_file("vm.ses")')

def initialize():
    h.finitialize()
    h.run()
    
initialize()

#save files
np.savetxt('01_vm_soma.txt', np.column_stack((np.array(cell.time_vector), np.array(cell.vm_soma))), delimiter = ' ')

img = plt.plot(np.array(cell.time_vector), np.array(cell.vm_soma))
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.savefig('01_vm_soma.eps')
plt.close()

quit()
