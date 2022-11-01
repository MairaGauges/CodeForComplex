from InputParser import parse_input
import os
from readVTK import initializeFromVTK, create_grid, get_boundary_data_dim2, inlet_cells
import numpy as np
import pyvista as pv
import Ofpp_patched
from OFScraper import get_boundary_data
from Initializer  import initialize
from EddyDetection import detect_eddy
import pyvista as pv


case_name = 'B4oxy30'



inlets, outlets, y_lim, z_lim, angle, pressure, radial_dir, axial_dir, T_threshold, dim2 = parse_input(os.path.join(f"{os.getcwd()}", f"data", f"{case_name}", f"CRNB_input.dic"))


Ny, Nz, y, z, V, vx, vy, vz, T, rho, filterArray,interpolated, grid, clipped, sl, VTK, domain = initializeFromVTK(case_name,y_lim,z_lim,radial_dir,axial_dir)





eddy_id, eddy_list =  detect_eddy(case_name, Ny, Nz, y, z, vy, vz, filterArray)


print(eddy_list)
np.save('filter_Array_coarse.npy', filterArray)

np.save('eddy_id_coarse.npy', eddy_id)




