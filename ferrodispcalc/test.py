from PolarDisp import PolarLMP
from TypeMap import UniPero
from PolarPolter import DispPloter

#a = PolarLMP('../traj.lammpstrj', UniPero)
#a.summary()
#a.get_polar('test',ele=['Ti'],nn_ele=['O'],r=3.7,num=6,vacancy=True,save=True)

b = DispPloter('test_disp.npy','test_origin.npy')

b.plot_polarization_time('xy',3,[10,10,10],save=True,dir='test',prefix='test')