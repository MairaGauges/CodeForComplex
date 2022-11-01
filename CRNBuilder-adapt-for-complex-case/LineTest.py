import numpy as np 
import pyvista as pv
from InputParser import parse_input
import os
from readVTK import initializeFromVTK, create_grid, get_boundary_data_dim2, inlet_cells



case_name='B4oxy30'

inlets, outlets, y_lim, z_lim, angle, pressure, radial_dir, axial_dir, T_threshold, dim2 = parse_input(os.path.join(f"{os.getcwd()}", f"data", f"{case_name}", f"CRNB_input.dic"))


grid = create_grid(z_lim,y_lim,radial_dir, axial_dir)[0]
line1 = pv.Line([0, 0.014, 0.0534999999999465],[0, 0.017, 0.0534999999999465])
line2 = pv.Line([0, 0.0225, 0.0534999999999465], [0, 0.024, 0.0534999999999465])
line4 = pv.Line([0, 0.0,-0.5999999], [0,0,-0.2])

list = grid.find_cells_along_line([0, 0.014, 0.0534999999999465], [0, 0.017, 0.0534999999999465])


print(list)


list2 = grid.find_cells_along_line([0, 0.0225, 0.0534999999999465], [0, 0.024, 0.0534999999999465])
print(list2)


list3= grid.find_cells_along_line([0, 0.19, 0.024999999999975], [0, 0.2, 0.024999999999975])
print(list3)


list4= grid.find_cells_along_line([0, 0.0,-0.5999999], [0,0,-0.2])
print(list4)


p = pv.Plotter()

p.add_mesh(grid, show_edges = True, opacity = 0.3)


p.add_mesh(line1, color = 'red')
p.add_mesh(line2, color = 'red')
p.add_mesh(line4, color = 'green')

#p.show_grid()
p.show()







