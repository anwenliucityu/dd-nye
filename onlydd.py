from copy import deepcopy
import datetime
import numpy as np
import matplotlib.pyplot as plt
import atomman as am
import atomman.unitconvert as uc


#base_system = am.load('atom_data','files/simple_cubic_edge_for_dd_perfect_reference.dat')
#disl_system = am.load('atom_data', 'files/simple_cubic_edge_for_dd.dat')
base_system = am.load('atom_data','files/screw_10_10_3/simple_cubic_screw_for_dd_perfect_reference.dat')
disl_system = am.load('atom_data', 'files/screw_10_10_3/simple_cubic_screw_for_dd.dat')
#base_system = am.load('atom_data', 'files/Al_edge_nye_per_ref.dat')
#disl_system = am.load('atom_data', 'files/Al_edge_nye_dislo.dat')
#base_system = am.load('atom_dump', 'files/fcc_Al_base.dump')
#disl_system = am.load('atom_dump', 'files/fcc_Al_disl.dump')
#base_system = am.load('atom_data', 'files/nye_per_ref.dat')
#disl_system = am.load('atom_data', 'files/nye_dislo.dat')

base_system.pbc=(True,True,True)
disl_system.pbc=(True,True,True)

alat=4.0
#alat = 4.05000466178543
cut_index=1.3
neighbors = base_system.neighborlist(cutoff = cut_index*alat)
dd = am.defect.DifferentialDisplacement(base_system, disl_system, neighbors=neighbors, reference=0)

ddmax =4.05/2**0.5/4
ddmax =2
print(ddmax)
params = {}
params['plotxaxis'] = 'x'
params['plotyaxis'] = 'y'
#params['xlim'] = (-12,-2)
params['xlim'] = (-12,12)
params['ylim'] = (-12,12)
#params['zlim'] = (-0.01, alat*6**0.5 / 2 + 0.01) # Should be one periodic width (plus a small cushion)
#params['zlim'] = (-0.1,3)
#params['zlim'] = (-0.01, 1.44)
params['figsize'] = 10         # Only one value as the other is chosen to make plots "regular"
params['arrowwidth'] = 1/50    # Made bigger to make arrows easier to see
params['arrowscale'] = 1.5#2.4     # Typically chosen to make arrows of length ddmax touch the corresponding atom circles
params['atomcmap'] = 'gray'
params['figsize'] = 20
#params['alat'] = 0 #to dismiss some ddvector if less than 0.1*burers

dd.ddplot('y', ddmax, **params)
plt.title('DD[2] component', fontsize=40)


# Example #1: Give base and disl systems with a cutoff
strain = am.defect.Strain(base_system, basesystem=disl_system, cutoff=cut_index*alat, theta_max=27)
#strain = am.defect.Strain(disl_system, basesystem=base_system, cutoff=1.25*    alat)                                                                           
strain.save_to_system()
straindict = strain.asdict()
#### plot seetings ##
par={}
par['fontsize'] =40
#####
####params for nye###
params['cmap'] = 'bwr'
params['xbins'] = 200
params['ybins'] = 200
params['scale'] = 0.1/0.09                         
params['nye_index'] = [2,2]

nye = 'Nye'+'['+str(params['nye_index'][0]+1)+','+str(params['nye_index'][1]+1)+']'
params['alat'] = 1000 #to dismiss some ddvector if less than 0.1*burers
dd.plot(base_system, 'x', ddmax, **params)
plt.title(nye,**par)
plt.show()  

