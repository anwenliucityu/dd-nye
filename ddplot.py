from copy import deepcopy
import datetime
import numpy as np
import matplotlib.pyplot as plt
import atomman as am
import atomman.unitconvert as uc

base_system = am.load('atom_data','/gauss12/home/cityu/anwenliu/git/dd/XMEAM/nye_per_ref.dat')
disl_system = am.load('atom_data','/gauss12/home/cityu/anwenliu/git/dd/XMEAM/nye_dislo.dat')

base_system.pbc=(False,False,True)
disl_system.pbc=(False,False,True)


#alat = uc.set_in_units(4.0, 'Å')
#alat = uc.set_in_units(2.93749682368896, 'Å') Ti
alat = uc.set_in_units(3.2056898354662, 'Å') 

cut_index = 1.3
neighbors = base_system.neighborlist(cutoff = cut_index*alat)
dd = am.defect.DifferentialDisplacement(base_system, disl_system, neighbors=neighbors, reference=1) #dislocated frame
# dd = am.defect.DifferentialDisplacement(base_system, disl_system, neighbors=neighbors, reference=0) #perfect frame

#ddmax = 2.94749682368896/2   # a/2<110> fcc dislocations use |b|/4 4 for partial, 2 for full
ddmax = np.sqrt((3.2056898354662)**2 + 5.20332683029676**2)/4
#ddmax = np.sqrt(4.65842**2 + 2.94749682368896**2)/2

# Set dict of keyword parameter values (just to make settings same for all plots below)
params = {}
params['plotxaxis'] = 'x'
params['plotyaxis'] = 'y'
params['xlim'] = (-10,25)
params['ylim'] = (-10,10)

#params['zlim'] = (-0.01, alat*6**0.5 / 2 + 0.01) # Should be one periodic width (plus a small cushion)
#params['zlim'] = (-0.1,1.5)
#params['zlim'] = (-0.01, 1.44)
params['figsize'] = 10         # Only one value as the other is chosen to make plots "regular"
params['arrowwidth'] = 1/50    # Made bigger to make arrows easier to see
params['arrowscale'] = 1.2#2.4     # Typically chosen to make arrows of length ddmax touch the corresponding atom circles
params['atomcmap'] = 'gray'
params['figsize'] = 20
params['alat'] = alat #to dismiss some ddvector if less than 0.1*burers

####params for nye###
params['cmap'] = 'bwr'
params['xbins'] = 200
params['ybins'] = 200
params['scale'] = 1
params['nye_index'] = [2,2]

#strain = am.defect.Strain(base_system, basesystem=disl_system, cutoff=cut_index*alat, theta_max=55)
strain = am.defect.Strain(disl_system, basesystem=base_system, cutoff=cut_index*alat, theta_max=55)
#strain = am.defect.Strain(disl_system, basesystem=base_system, cutoff=1.25*alat) 
strain.save_to_system()
straindict = strain.asdict()

#### plot seetings ##
par={}
par['fontsize'] = 40
#####

nye = 'Nye'+'['+str(params['nye_index'][0]+1)+','+str(params['nye_index'][1]+1)+']'
fig = dd.plot(disl_system, 'z', ddmax, **params)
#print(fig.x0, fig.y0, fig.x_sq, fig.y_sq, fig.xy)
plt.title(nye+' with DD[3]',**par)
plt.xlabel('x ('+ u'\u212B'+')', fontsize = 45)
plt.ylabel('y ('+ u'\u212B'+')', fontsize = 45).set_rotation(0)
plt.show()
