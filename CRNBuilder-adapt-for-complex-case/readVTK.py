
import myUtils as mU
import numpy as np
import pyvista as pv
import math as m
import os

def get_boundary_data_dim2(inlets, radial_dir, axial_dir, y_lim, z_lim, T, vx, vy, vz, rho, y, z, Ny, Nz, filterArray):

    '''Returns the dictionary boundary data which contains a dictionary for every inlet. 
       Each inlet dictionary contains arrays for the properties U, T and rho at the inlet. The information is extracted from the cell information.'''

    #create dirctionary boundary_data
    boundary_data = {}
    properties = ['U','T', 'rho']
    inlet_cells_ = inlet_cells(inlets, y_lim, z_lim, radial_dir, axial_dir)
    print(inlet_cells_)
    #create empty arrays--> need number of cells
    for boundary_name in inlet_cells_ :
        boundary_data[boundary_name] = {}
        cells = []
        cell_list = inlet_cells_[boundary_name]
        num_cells = len(cell_list)
        

        for k in range(num_cells):

           j = int (cell_list[k] % Nz)
           i = int ((cell_list[k] - j)/Nz)
 
           
           #j  = (cell_list[k] % Ny) -1
           #i = int((cell_list[k] - (j+1))/Ny)
           print(i,j)
           if filterArray[i,j] == 1:
               cells.append((i,j))

        array = 0
        for prop in properties:
            if prop == 'T':
                array = T
            elif prop == 'rho':
                array = rho
           
            if prop == 'T' or prop == 'rho':
                new_array = np.empty(shape= (len(cells)))
                for k in range(len(cells)):
                   i,j = cells[k]
                   new_array[k] = array[i,j]
            elif prop == 'U':
                new_array = np.empty(shape = (len(cells),3))
                for k in range(len(cells)):
                    i,j = cells[k] 
                    new_array[k,0] = vx[i,j]
                    new_array[k,1] = vy[i,j]
                    new_array[k,2] = vz[i,j]

            boundary_data[boundary_name][prop] = new_array
             
    return boundary_data
       


 
def inlet_cells(inlets, y_lim, z_lim, radial_dir, axial_dir):
    '''Returns the dictionary inlet_cells. This dictionary has a key for every inlet. The corresponding value is a list of the cells indecis (from stuctured grid) that are part of the grid. '''

    grid = create_grid(z_lim,y_lim, radial_dir, axial_dir)[0]
    inlet_cells = {}
    
    point1 = [0,0,0]
    point2 = [0,0,0]
    for i in range(len(inlets)):
        inlet_name = inlets[i][0]
        direction = inlets[i][1][0] 
        start = float(inlets[i][1][1]) 
        end = float(inlets[i][1][2]) 
        loc = float(inlets[i][1][3])*(1-(1e-12))
        #radial        
        if direction == 0:
            point1[radial_dir] = start
            point2[radial_dir] = end 
            point1[axial_dir] = loc
            point2[axial_dir] = loc
        
        #axial
        elif direction == 1:
            point1[radial_dir] = loc
            point2[radial_dir] = loc
            point1[axial_dir] = start
            point2[axial_dir] = end
       
        print((point1,point2)) 
        cell_list = grid.find_cells_along_line(point1,point2)
        print(cell_list)
        inlet_cells[inlet_name] = cell_list
    return inlet_cells


def create_grid(z_lim,y_lim,radial_dir, axial_dir):
    '''A grid is created according to the size of the original domain (given in input file). The direction of the grid depends on the what the radial and axial direction of the original domain are.'''

    #create new mesh
    zmin = z_lim[0]
    zmax = z_lim[1]
    ymin = y_lim[0]
    ymax = y_lim[1]

    #decide how to define an appropriate cell size
    #B4oxy30 case
    #stepY = 0.001
    #stepZ = 0.001307
    stepY = 0.002
    stepZ = 0.005228   



    #HM1 bluff body case 
    #stepY = 0.00015
    #stepZ = 0.0012    

    z  =np.arange(zmin,zmax+stepZ , stepZ)
    y = np.arange(ymin,ymax+stepY, stepY)

    y_t, z_t = np.meshgrid(y,z)
    #z_t, y_t = np.meshgrid(z, y) 
    x_t = np.zeros(z_t.shape)


    #Variable assignment needs to be adapted if not y radial and z axial
    # x --> None
    # y --> radial direction
    # z --> axial direction

    if axial_dir == 1:
        y = z_t
        if radial_dir == 0:
            x = y_t
            z = x_t
        elif radial_dir == 2:
            z = y_t
            x = x_t

    elif axial_dir == 0:
        x = z_t
        if radial_dir == 2:
            z = y_t
            y = x_t
        elif radial_dir == 1:
            z = x_t
            y = y_t

    elif axial_dir == 2:
        z = z_t
        if radial_dir == 0:
            x = y_t
            y = x_t
        elif radial_dir ==1:
            y = y_t
            x = x_t
           
    grid = pv.StructuredGrid(x,y,z)

    #find Ny and Nz from grid created
    Nz =int((zmax-zmin)/stepZ)
    Ny = int((ymax - ymin)/stepY)

    return grid,Nz,Ny



def initializeFromVTK(case_name,y_lim,z_lim, radial_dir, axial_dir):
    ''' Reads VTK data and slices and halves domain so that a 2d 'wedge' is left. Structured grid created and for every cell center of the structured grid it is checked if this point is part of the original domain. From this a filterArray is created indicating which cell is part of original domain and which is not. Properties are interpolated from original domain onto strutured grid and the properties data restructured into an array in teh shape of the structured grid. Cells that were not part of original domain are marked with NaN values.  '''

    #get path 
    path  = os.path.join(f"{os.getcwd()}",f"data/")
    folderName = case_name
  
    #read VTK data
    timeL,fileL = mU.readVTKSeries(path,folderName)
    #VTK is a pyvista pointset Unstructured grid
    VTK = pv.read(fileL[-1])
    
    directions = [radial_dir, axial_dir]
    dir_num = [0,1,2]
    dir_strings = ['x','y','z']
     
    sl_dir =dir_strings[(np.setdiff1d(dir_num,directions))[0]]
    clip_dir = '-'+dir_strings[radial_dir]    

    #clip domain 
    clipped = VTK.clip(clip_dir,invert = True)
    #slice is Polydata
    sl = VTK.slice(normal=sl_dir)
    
    #HM1 case
    #domain = VTK

    #complex case
    domain = sl.clip(clip_dir,invert = True)

    grid, Nz, Ny = create_grid(z_lim,y_lim,radial_dir, axial_dir)

    filterArray = np.empty(shape = (Ny,Nz))
 
    #get centerpoint of grid's cells
    cell_centers = grid.cell_centers().points

    #interpolate from VTK to structured grid - interpolated is a structured grid, point set

    #for comlex case
    interpolated = grid.interpolate(clipped, radius = 0.01, sharpness = 10)

    #for HM1 case
    #interpolated = grid.interpolate(sl, radius = 0.001, sharpness = 10)

    interpolated = interpolated.point_data_to_cell_data()

    #get arrays
    #switched from ...Mean
    pT = interpolated.get_array('TMean') 
    pU = interpolated.get_array('UMean')
    prho  = interpolated.get_array('rhoMean')
     
            

    pvy = []
    pvz = []
    pvx = []

    for k in range(len(pU)):  
        pvx.append(pU[k][0])
        pvy.append(pU[k][1])
        pvz.append(pU[k][2])

    if axial_dir == 0:
        vz_temp = pvx.copy()
    elif axial_dir == 1:
        vz_temp = pvy.copy()
    elif axial_dir == 2:
        vz_temp = pvz.copy()
    
    if radial_dir == 0:
        vy_temp = pvx.copy()
    elif radial_dir == 1:
        vy_temp = pvy.copy()
    elif radial_dir == 2:
        vy_temp = pvz.copy()

    #which direction is not axial and not radial
    Ndir = np.setdiff1d(dir_num,directions)[0]
    if Ndir == 0:
        vx_temp = pvx.copy()
    elif Ndir == 1:
        vx_temp = pvy.copy()
    elif Ndir == 2:    
        vx_temp = pvz.copy()
    #create empty arrays 
    V = np.empty(shape=(Ny,Nz))
    y = np.empty(shape=(Ny,Nz))
    z = np.empty(shape=(Ny,Nz))
    T = np.empty(shape=(Ny,Nz))
    rho = np.empty(shape=(Ny,Nz))
    vx = np.empty(shape=(Ny,Nz))
    vy = np.empty(shape=(Ny,Nz))
    vz  = np.empty(shape=(Ny,Nz))
    #restructuring data --> switch i and j since the build up array is opposite to the counting of cells with pyvista
    

    for i in range(Ny):
        for j in range(Nz):
            counter = Nz*i + j
            percent = round((counter/(Ny*Nz))*100, 2)
            print('restructuring of data completed to ' + str(percent)+ ' %', end = '\r')
            k = domain.find_containing_cell(cell_centers[counter])
            if k == -1:
                filterArray[i,j] = 0
                T[i,j] = None
                vx[i,j] = None
                vy[i,j] = None
                vz[i,j] = None
                rho[i,j] = None
                y[i,j] = None
                z[i,j] = None
                V[i,j] = None

            else:
                filterArray[i,j] = 1
                T[i,j] = pT[counter]
                vx[i,j] = vx_temp[counter]
                vy[i,j] = vy_temp[counter]
                vz[i,j] = vz_temp[counter]
                rho[i,j] = prho[counter]

                #assignments depend on direction
                y[i,j] = cell_centers[counter][radial_dir]
                z[i,j] = cell_centers[counter][axial_dir]
                bounds = interpolated.cell_bounds(counter)
                ri = bounds[2*radial_dir]
                ro = bounds[2*radial_dir+1]
                h = bounds[axial_dir*2]-bounds[axial_dir*2 +1]
                V[i,j] = m.pi*h*(ro**2 - ri**2)
    

    ''' Wrong!
    for j in range(Nz):
        for i in range(Ny):
            counter = Ny*j +i
            percent = round((counter/(Ny*Nz))*100, 2)
            print('restructuring of data completed to ' + str(percent)+ '%', end = '\r')
            k = domain.find_containing_cell(cell_centers[counter])
            #if k == -1:
            if k == 6:
                filterArray[i,j] = 0
                T[i,j] = None
                vx[i,j] = None
                vy[i,j] = None
                vz[i,j] = None
                rho[i,j] = None 
                y[i,j] = None
                z[i,j] = None
                V[i,j] = None

            else:
                filterArray[i,j] = 1
                T[i,j] = pT[counter]
                vx[i,j] = vx_temp[counter]
                vy[i,j] = vy_temp[counter]
                vz[i,j] = vz_temp[counter]
                rho[i,j] = prho[counter] 
              
                #assignments depend on direction
                y[i,j] = cell_centers[counter][radial_dir]
                z[i,j] = cell_centers[counter][axial_dir]
                bounds = interpolated.cell_bounds(counter)
                #x x y y z z 
                ri = bounds[2*radial_dir]
                ro = bounds[2*radial_dir+1]
                h = bounds[axial_dir*2]-bounds[axial_dir*2 +1]
                V[i,j] = m.pi*h*(ro**2 - ri**2)
      '''

 
    return Ny, Nz, y, z, V, vx, vy, vz, T, rho, filterArray




