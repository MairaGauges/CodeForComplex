
from InputParser import parse_input
import os
from readVTK import initializeFromVTK, create_grid, get_boundary_data_dim2, inlet_cells
import numpy as np
import pyvista as pv
import Ofpp_patched
from OFScraper import get_boundary_data
#from CRNBuilder import post_process_boundary_data
from Initializer  import initialize 

from EddyDetection import detect_eddy
import pyvista as pv

case_name='B4oxy30'
#case_name = 'HM1_bluff-body_flame'



inlets, outlets, y_lim, z_lim, angle, pressure, radial_dir, axial_dir, T_threshold, dim2 = parse_input(os.path.join(f"{os.getcwd()}", f"data", f"{case_name}", f"CRNB_input.dic"))


Ny, Nz, y, z, V, vx, vy, vz, T, rho, filterArray,interpolated, grid, clipped, sl, VTK, domain = initializeFromVTK(case_name,y_lim,z_lim,radial_dir,axial_dir)


#complex case
#Ny = 210
#Nz = 500

#HM1 case
#Ny = 899
#Nz = 275


'''
np.save('y.npy',y)
np.save('z.npy', z)
np.save('filterArray.npy',filterArray)
np.save('vz.npy',vz)
np.save('vy.npy',vy)
np.save('T.npy', T)
np.save('y.npy', y)
np.save('z.npy',z)
'''


'''
y = np.load('yH.npy')
z = np.load('zH.npy')
filterArray = np.load('filterArrayH.npy')
vz = np.load('vzH.npy')
vy = np.load('vyH.npy')
'''


filterArray = np.load('filterArray.npy')
y = np.load('y.npy')
z  = np.load('z.npy')
vz = np.load('vz.npy')
vy = np.load('vy.npy')

eddy_centers = [[32, 8], [7, 8], [186, 76], [162, 276], [141, 426], [33, 443]]



eddy_id, eddy_list =  detect_eddy(case_name, Ny, Nz, y, z, vy, vz, filterArray)

#print(eddy_list)
#eddy_id = detect_eddy(case_name, Ny, Nz, y, z, vy, vz)



np.save('eddy_id_big.npy',eddy_id)


'''
p = pv.Plotter()
p.add_mesh(domain, opacity = 0.4, color = 'green')
p.add_mesh(grid, color = 'red')
#p.add_mesh(sl, scalars = 'T')
#p.add_mesh(clipped, scalars = 'T')
#p.add_mesh(interpolated, scalars = 'T', cmap = 'jet')
p.show_grid()
p.show()
'''



'''
properties = ['U','T', 'rho']

for inlet in boundary_data:
    print(boundary_data[inlet])
    
    for p in properties: 
        print(p)
        print(boundary_data[inlet][p])

'''



'''
p = pv.Plotter()
p.add_mesh(domain)
p.show_grid()
p.show()

'''


