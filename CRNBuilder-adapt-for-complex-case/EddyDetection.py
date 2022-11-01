from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
import os

from Initializer import initialize


def detect_eddy(case_name, Ny, Nz, y, z, vy, vz, filterArray,i_min=None, j_min=None, i_max=None, j_max=None, a=4, b=3, extended_output=False, load_cached=True, save_cache=True):
    print("Detecting eddies...")


    if load_cached and not extended_output:
        path = os.path.join(f"{os.getcwd()}", "data", f"{case_name}", "cache", "eddy")
        if os.path.isdir(path):
            print("Cache found. Loading data...")
            eddies_ij = []
            for file in os.listdir(path):
                eddy_i = []
                eddy_j = []
                with open(os.path.join(f"{path}", f"{file}"), 'r') as f:
                    lines = f.readlines()
                for line in lines:
                    eddy_i.append(int(line.split(',')[0].strip(' ')))
                    eddy_j.append(int(line.split(',')[1].rstrip('\n').strip(' ')))
                eddies_ij.append(zip(eddy_i, eddy_j))
            eddy_id = np.zeros((Ny, Nz), dtype=int)
            for e in range(len(eddies_ij)):
                for i, j in eddies_ij[e]:
                    eddy_id[i, j] = e + 1
            return eddy_id
        else:
            print("Cache not found. Searching eddies...")

    if i_min is None:
        i_min = a + 1
    if j_min is None:
        j_min = a + 1
    if i_max is None:
        i_max = y.shape[0] - a - 1
    if j_max is None:
        j_max = y.shape[1] - a - 1

    print("Applying 1st constraint...")
    # First constraint
    i_o = []
    j_o = []
    for j in range(j_min, j_max):
        for i in range(i_min, i_max - 1):
            #check in square area around it if any cells are outside of orginal domain --> if so do not consider as eddy center, potential eddy centrs must be far enough from border
           
            outside = False
            for k in range(i-a,i+a+1):
                for m in range(j-a,j+a):
                    if filterArray[k,m] == 0:
                        outside = True
            if outside:
                continue

            if vz[i, j] * vz[i + 1, j] < 0 and \
                    vz[i - a, j] * vz[i + 1 + a, j] < 0 and \
                    abs(vz[i, j]) < abs(vz[i - a, j]) and \
                    abs(vz[i + 1, j]) < abs(vz[i + a, j]):
                i_o.append([i, i + 1])
                j_o.append([j, j])

    # Second constraint
    print("Applying 2nd constraint...")
    new_i_o = []
    new_j_o = []
    for i_pair, j_pair in zip(i_o, j_o):
        # TODO: Add check for vy sign (cf. https://doi.org/10.1175/2009JTECHO725.1 pg.5 second constraint)
        #  Everything seems to work even without the check, so I'll leave it there, however if there are problems in the eddy detection phase this fix should be one of the first thing to try

        if abs(vy[i_pair[0], j_pair[0]]) < abs(vy[i_pair[0], j_pair[0] + a]) and \
                abs(vy[i_pair[1], j_pair[1]]) < abs(vy[i_pair[1], j_pair[1] + a]) and \
                abs(vy[i_pair[0], j_pair[0]]) < abs(vy[i_pair[0], j_pair[0] - a]) and \
                abs(vy[i_pair[1], j_pair[1]]) < abs(vy[i_pair[1], j_pair[1] - a]):
            new_i_o.append(i_pair)
            new_j_o.append(j_pair)

    i_o = [k for row in new_i_o for k in row]
    j_o = [k for row in new_j_o for k in row]

    # Third constraint (velocity local minimum)
    # TODO: There's 100% a smarter way to do this but it works so I'll leave it like that
    print("Applying 3rd constraint...")
    vmin_coords = []
    for i, j in zip(i_o, j_o):
        vmin_coords.append([i, j])
        for u in range(i - b, i + b + 1):
            for v in range(j - b, j + b + 1):
                #checking if new minimum is far enough away from border of original domain
                outside = False

                for k in range (u-a, u+a+a):
                    for m in range(v-a,v+a+1):
                        if k > Ny-1 or k < 0 or m < 0 or m > Nz-1 or filterArray[k,m] == 0:
                            outside = True
                if outside:
                    continue            

                if velocity_magnitude(u, v, vy, vz) < velocity_magnitude(vmin_coords[-1][0], vmin_coords[-1][1], vy, vz):
                    vmin_coords[-1] = [u, v]
        vmin_coords_temp = vmin_coords[-1]


        for it in range(10000):
            for u in range(vmin_coords[-1][0] - b, vmin_coords[-1][0] + b + 1):
                for v in range(vmin_coords[-1][1] - b, vmin_coords[-1][1] + b + 1):
                    outside = False

                    for k in range (u-a, u+a+a):
                        for m in range(v-a,v+a+1):
                            if k > Ny-1 or k < 0 or m < 0 or m > Nz-1 or filterArray[k,m] == 0:
                                outside = True

                    if outside:
                        continue   
                    
                    if velocity_magnitude(u, v, vy, vz) < velocity_magnitude(vmin_coords[-1][0], vmin_coords[-1][1], vy, vz):
                        vmin_coords[-1] = [u, v]
            if vmin_coords[-1] == vmin_coords_temp:
                break
    vmin_coords = [ii for n, ii in enumerate(vmin_coords) if ii not in vmin_coords[:n]]
    print(f"Found {len(vmin_coords)} possible eddies.")

    # Fourth constraint (the velocity vectors around the center rotate coherently)
    print("Applying 4th constraint...")
    to_keep = []
    for coords in vmin_coords:
        to_keep.append(True)
        quad = []
        i, j = coords
        # Starting from the southwest corner proceed counterclockwise
        v = j - (a - 1)
        for u in range(i - (a - 1), i + (a - 1)):
            quad.append(get_quadrant(vy[u, v], vz[u, v]))
        u = i + (a - 1)
        for v in range(j - (a - 1), j + (a - 1)):
            quad.append(get_quadrant(vy[u, v], vz[u, v]))
        v = j + (a - 1)
        for u in range(i + (a - 1), i - (a - 1), -1):
            quad.append(get_quadrant(vy[u, v], vz[u, v]))
        u = i - (a - 1)
        for v in range(j + (a - 1), j - (a - 1) - 1, -1):
            quad.append(get_quadrant(vy[u, v], vz[u, v]))
        for q in range(len(quad) - 1):
            if quad[q] % 4 == quad[q + 1] % 4 or quad[q] % 4 == (quad[q + 1] - 1) % 4:
                continue
            to_keep[-1] = False

    vmin_coords = [vmin_coords[ii] for ii in range(len(vmin_coords)) if to_keep[ii]]
    print(f"Found {len(vmin_coords)} eddies.")

    #vmin_coords = [[32, 8], [7, 8], [186, 76], [162, 276], [141, 426], [33, 443]]
    #vmin_coords = [[141,426]]


    # Get actual entire eddy
    paths_y = []
    paths_z = []
    paths_ij = []
    eddies_ij = []


    

    for coords in vmin_coords:
        print(f"Searching the contour of the eddy number {vmin_coords.index(coords)}...")
        i, j = coords
        closed = [True, True]
        k = 1
        paths_y_temp = []
        paths_z_temp = []
        paths_i_temp = []
        paths_j_temp = []
        stopped = False
        #davor bis -11 bzw pop bei 12
        while closed[-1] or closed[-2] or closed[-3] or closed[-4] or closed[-5] or closed[-6] or closed[-7] or closed[-8] or closed[-9] or closed[-10] or closed[-11] or closed[-12] or closed[-13] or closed[-14] or closed[-15]:
       
        #while closed[-1] or closed[-2] or closed[-3] or closed[-4] or closed[-5]:


            
            #if outside of domain do not continue looking for eddy (does this make eddy smaller?)
            if (j+k)>= Nz or filterArray[i,j+k] == 0: 
                count = 0
                print('stopped')
                print((closed)) 
                stopped = True
                for cl in reversed(closed):
                    if cl: 
                        break
                    count += 1  
                #counter back to 0 if closed loop inbetween ...               
                num_true = len(paths_y_temp) - count -1
                for m in range(15-count):
                    closed.append(False)
               
                print('num_true = ' + str(num_true))

                print(len(paths_y_temp)) 
                
                paths_y_temp = paths_y_temp[num_true:]
                paths_i_temp = paths_i_temp[num_true:]
                paths_z_temp = paths_z_temp[num_true:]
                paths_j_temp = paths_j_temp[num_true:]
                #paths_j_temp = paths_j_temp[num_true-3:] 
                print(len(paths_y_temp))
 
                break 



            path_y, path_z, path_i, path_j, close = simulate_particle(i, j + k, y, z, vy, vz,filterArray)
            closed.append(close)


            paths_y_temp.append(path_y)
            paths_i_temp.append(path_i)
            if len(paths_y_temp) > 16:
                paths_y_temp.pop(0)
                paths_i_temp.pop(0)
            paths_z_temp.append(path_z)
            paths_j_temp.append(path_j)
            if len(paths_z_temp) > 16:
                paths_z_temp.pop(0)
                paths_j_temp.pop(0)
            k += 1


            if stopped:
                break

        print(closed)
        paths_y.append([paths_y_temp[-1]])
        paths_y.append([paths_y_temp[-1]])





        #paths_z.append([paths_z_temp[0], paths_z_temp[2]])
        #paths_y.append([paths_y_temp[0], paths_y_temp[2]])


        paths_ij.append(zip(paths_i_temp[0], paths_j_temp[0]))
        

        print(f"Filling the eddy number {vmin_coords.index(coords)}...")
        eddies_ij.append(find_enclosed_points(paths_i_temp[0], paths_j_temp[0]))
        

    # Get coordinates from i and j
    #y_o = [y[i, j] for i, j in zip(i_o, j_o)]
    #z_o = [z[i, j] for i, j in zip(i_o, j_o)]
    #y_vmin = [y[vmin_coords[i][0], vmin_coords[i][1]] for i in range(len(vmin_coords))]
    #z_vmin = [z[vmin_coords[i][0], vmin_coords[i][1]] for i in range(len(vmin_coords))]
    
    '''if save_cache:
        print("Saving cache...")
        path = os.path.join(f"{os.getcwd()}", f"data", f"{case_name}", f"cache", f"eddy")
        if not os.path.isdir(path):
            os.mkdir(path)
        else:
            for file in os.listdir(path):
                os.remove(file)
        for e in range(len(eddies_ij)):
            to_file = ""
            for i, j in eddies_ij[e]:
                to_file += f"{i}, {j}\n"
            with open(os.path.join(f"{path}", f"{e}.csv"), 'w') as f:
                f.write(to_file)'''

    if not extended_output:
        eddy_id = np.zeros((Ny, Nz), dtype=int)
        for e in range(len(eddies_ij)):
            for i, j in eddies_ij[e]:
                eddy_id[i, j] = e + 1
        return eddy_id
    #else:
     #   return y_o, z_o, y_vmin, z_vmin, paths_y, paths_z, paths_ij, eddies_ij


def velocity_magnitude(i, j, vy, vz):
    return sqrt(vy[i, j] * vy[i, j] + vz[i, j] * vz[i, j])


def get_quadrant(v1, v2):
    if v1 >= 0 and v2 >= 0:
        return 0
    elif v1 < 0 and v2 >= 0:
        return 1
    elif v1 < 0 and v2 < 0:
        return 2
    elif v1 >= 0 and v2 < 0:
        return 3
    raise RuntimeError("The quadrant can not be identified.")


def simulate_particle(i0, j0, y, z, vy, vz, filterArray,dt=1e-5, max_iter=int(1e5), precision=1e-3):
    # print(f"Simulating particle at {i0 = } and {j0 = }")
    #dt originally 1e-5
    # precision originally 1e-3
    # max iter originalyy 1e5
    y_min = np.nanmin(y)
    y_max = np.nanmax(y)
    z_min = np.nanmin(z)
    z_max = np.nanmax(z)
    par_y0 = y[i0, j0]
    par_z0 = z[i0, j0]
    par_y = par_y0
    par_z = par_z0
    path_y = [par_y]
    path_z = [par_z]
    i = i0
    j = j0
    path_i = [i]
    path_j = [j]
    close = False 
    close_1 = False
    close_2 = False
    while i == i0 and j == j0:
        time = dt
        outside_domain = True
        
        while outside_domain:
            par_y += vy[i, j] * time
            par_y = max(min(par_y, y_max), y_min)
            par_z += vz[i, j] * time
            par_z = max(min(par_z, z_max), z_min)
            i = np.argmin(np.abs(y[:, 0] - par_y))
            j = np.argmin(np.abs(z[0, :] - par_z))
            if filterArray[i,j] == 1:
                outside_domain = False
                path_y.append(par_y)
                path_z.append(par_z)
            time /=1.25            
        add_ij = True
        for k in range(len(path_i)):
            if path_i[k] == i and path_j[k] == j:
                add_ij = False
                break
        if add_ij:
            path_i.append(i)
            path_j.append(j)

    #condition needed if particle has moved far enough away from starting point 
    away = [False, False, False]
    for t in range(max_iter):
        time = dt
        outside_domain = True
 
        while outside_domain:
            par_y += vy[i, j] * time
            par_y = max(min(par_y, y_max), y_min)
            par_z += vz[i, j] * time
            par_z = max(min(par_z, z_max), z_min)
            i = np.argmin(np.abs(y[:, 0] - par_y))
            j = np.argmin(np.abs(z[0, :] - par_z))
            if filterArray[i,j] == 1:
                outside_domain = False
                path_y.append(par_y)
                path_z.append(par_z)
            time /=1.25        
        add_ij = True
        for k in range(len(path_i)):
            if path_i[k] == i and path_j[k] == j:
                add_ij = False
                break
        if add_ij:
            path_i.append(i)
            path_j.append(j)
        #try other condition for closed path?? --> look at it cellwise, are i and j close enough to i0 and j0?


        '''
        if abs(par_y - par_y0) <= precision * (y_max - y_min) and \
                abs(par_z - par_z0) <= precision * (z_max - z_min):
            close = True
            break
        '''
        #'''
        if not i in range(i0-1,i0+2) or not j in range(j0-2,j0+3):
            away[0] = True

        if not i in range(i0-2,i0+3) or not j in range(j0-2,j0+3):
            away[1] = True

        if not i in range(i0-3,i0+4) or not j in range(j0-3,j0+4):
            away[2] = True

        if i in range(i0-1,i0+2) and j in range(j0-1,j0+2) and away[0]:
            close = True
            break 
        
        if i in range(i0-2,i0+3) and j in range(j0-2,j0+3) and away[1] and not close_1:
            close_1 = True
            end_point_1 = (i,j)
            index_1 = len(path_i)
        if i in range(i0-3,i0+4) and j in range(j0-3,j0+4) and away[2] and not close_2:
            close_2 = True
            end_point_2 = (i,j)
            index_2 = len(path_i)


    #close if hole of 1 cell --> fill these cells
    #currently connecting to starting point - maybe connect to closest point?
    if not close: 
        if close_1:
            close = True
            end_i,end_j = end_point_1[0],end_point_1[1]
            i = int(i0+0.5*(end_i-i0)) 
            j = int(j0 + 0.5*(end_j-j0))
            path_i = path_i[:index_1]
            path_j = path_j[:index_1]
            path_i.append(i)
            path_j.append(j)


        '''
        elif close_2: 
            close = True
            end_i,end_j = end_point_2[0],end_point_2[1]
            i = i0+(end_i-i0)/3
            j = j0 + (end_j-j0)/3
            path_i = path_i[:index_2]
            path_j = path_j[:index_2]
            path_i.append(i)
            path_j.append(j)
            i = i0+(end_i-i0)*2/3
            j = j0 + (end_j-j0)*2/3
            path_i.append(i)
            path_j.append(j)
        '''
        # if sum(abs(np.array(path_y[:-3]) - par_y) <= precision * (np.max(y) - np.min(y))) > 5 and \
        #         sum(abs(np.array(path_z[:-3]) - par_z) <= precision * (np.max(z) - np.min(z))) > 5:
        #     break
    return path_y, path_z, path_i, path_j, close


def find_enclosed_points(path_i, path_j):
    i_min = min(path_i)
    i_max = max(path_i)
    j_min = min(path_j)
    j_max = max(path_j)
    path_ij = []
    for i, j in zip(path_i, path_j):
        path_ij.append([i, j])
    outer = np.zeros((i_max + 1 - i_min, j_max + 1 - j_min), dtype=bool)
    eddy_i = []
    eddy_j = []
    # North and south boundaries
    for i in range(i_min, i_max + 1):
        if not [i, j_min] in path_ij:
            outer[i - i_min, j_min - j_min] = True
        if not [i, j_max] in path_ij:
            outer[i - i_min, j_max - j_min] = True
    # East and west boundaries
    for j in range(j_min, j_max + 1):
        if not [i_min, j] in path_ij:
            outer[i_min - i_min, j - j_min] = True
        if not [i_max, j] in path_ij:
            outer[i_max - i_min, j - j_min] = True
    # Central points
    for i in range(np.prod(outer.shape)):
        outer_new = outer.copy()
        for i in range(i_min + 1, i_max):
            for j in range(j_min + 1, j_max):
                if [i, j] in path_ij:
                    continue
                if outer[i - 1 - i_min, j - j_min] or \
                        outer[i + 1 - i_min, j - j_min] or \
                        outer[i - i_min, j - 1 - j_min] or \
                        outer[i - i_min, j + 1 - j_min]:
                    outer_new[i - i_min, j - j_min] = True
        if (outer_new == outer).all():
            break
        outer = outer_new
    # Build eddy_ij
    for i in range(outer.shape[0]):
        for j in range(outer.shape[1]):
            if not outer[i, j]:
                eddy_i.append(i + i_min)
                eddy_j.append(j + j_min)
    eddy_ij = zip(eddy_i, eddy_j)

    return eddy_ij


def main():
    # Load data
    case_name = "CH4Piloted"
    Ny, Nz, y, z, V, vx, vy, vz, T, rho = initialize(case_name)

    a = 4
    b = 3

    y_o, z_o, y_vmin, z_vmin, paths_y, paths_z, paths_ij, eddies_ij = detect_eddy(case_name, Ny, Nz, y, z, vy, vz, a=a, b=b, extended_output=True, save_cache=False)

    # Plots
    print("Plotting...")
    fig = plt.figure()

    ax = fig.add_subplot(121)
    # ax.streamplot(y.T, z.T, vy.T, vz.T, density=2, linewidth=0.5, arrowsize=0.5, color='black')
    for p in range(len(paths_y)):
        ax.plot(paths_y[p], paths_z[p], color='blue')
    ax.scatter(y_vmin, z_vmin, marker='x', color='blue')
    ax.set_title('Velocity streamlines')
    ax.set_xlabel('y')
    ax.set_ylabel('z')

    fig = plt.figure()

    ax = fig.add_subplot(111)
    for i in range(vy.shape[0]):
        for j in range(vy.shape[1]):
            vy_0 = vy[i, j]
            vz_0 = vz[i, j]
            vy[i, j] /= (vy_0 + vz_0)
            vz[i, j] /= (vy_0 + vz_0)
    ax.quiver(y, z, vy, vz, units='xy', angles='xy', scale_units='xy', scale=0.01)
    for p in range(len(paths_y)):
        ax.plot(paths_y[p], paths_z[p], color='blue')
        for i, j in eddies_ij[p]:
            ax.scatter(y[i, j], z[i, j], facecolors='none', edgecolors=[0, 0, 1, 0.5])
        for i, j in paths_ij[p]:
            ax.scatter(y[i, j], z[i, j], facecolors='none', edgecolors='red')
    ax.scatter(y_vmin, z_vmin, marker='x', color='red')
    ax.scatter(y_o, z_o, facecolor='none', edgecolor='red')
    ax.set_title('Velocity vectors')
    ax.set_xlabel('y')
    ax.set_ylabel('z')

    plt.show()


if __name__ == '__main__':
    main()
