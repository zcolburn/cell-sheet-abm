# -*- coding: utf-8 -*-
"""
Created on Sun May 27 19:54:25 2018

@author: ZacharyColburn
"""

# Import modules.
import cell
from paramup_fxn import paramup_fxn

import numpy as np
import matplotlib.pyplot as plt
from itertools import chain
import csv
import os

'''
Set parameters.
'''
'''
Name of the csv file to output that shows wound closure progress. Each column
represents one unit width in the environment. Each row represents a different
timepoint. Each cell is filled with the number of cells in that column of the
wound environment.
'''
filename = 'rJEB.csv'

num_cells = 2000# Cells to seed the environment with.
num_timepoints = 40
frame_pic_mod = 1# The number of timepoints between exporting each wound image.

proportion = 1# = rJEB / (rJEB + JEB) # For altering the ratio of rJEB:JEB

scratch_time = 10# Timepoint for scratching

'''
Default cell type parameters.

rjeb = cells with beta-4 integrin
jeb = cells without beta-4 integrin

mov_chg = randomness factor for changing x and y displacement vectors
cell_rad = cell radius
max_speed = maximum cell speed
division_rate = probability that a cell divides at a given timepoint
death_rate = probability that a cell dies at a given timepoint
min_division_dist = In order for this cell to divide its distance to neighbors
    must exceed this value.
div_steps = Number of orientations to try dividing
'''
rjeb_params = {'b4_integrin':1,
             'paxillin':1,
             'mov_chg':0.25,
             'cell_rad':1,
             'max_speed':1,
             'division_rate':0.05,
             'death_rate':0.01,
             'min_division_dist':0.5,
             'div_steps':8}

jeb_params = {'b4_integrin':0,
              'paxillin':1,
             'mov_chg':0.25,
             'cell_rad':1,
             'max_speed':1,
             'division_rate':0.05,
             'death_rate':0.01,
             'min_division_dist':0.5,
             'div_steps':8}

# Dimensions of the grid environment
x_grid_dim = 40
y_grid_dim = 20

# The bounds of the scratch
scratch_left = 7
scratch_right = 14


'''
Main
'''
# Set random seed.
np.random.seed(8)

# Grid dimensions
extra = {'xmin':0,
         'xmax':x_grid_dim,
         'ymin':0,
         'ymax':y_grid_dim}

'''
Create a grid. Each element will be filled with a unique id.

grid_indices = a matrix filled with unique ids
substrate = a dictionary with keys equal to the unique ids in grid_indices.
    For storing local substrate data.
medium = a dictionary with keys equal to the unique ids in grid_indices.
    For storing local medium data.
grid = a dictionary with keys equal to the unique ids in grid_indices.
    For storing local cell ids.
'''
grid_indices = np.zeros((x_grid_dim + 1, y_grid_dim + 1))

i = 0
substrate = {}
for x in range(len(grid_indices)):
    for y in range(len(grid_indices[0])):
        grid_indices[x][y]  = i
        substrate[i] = {}
        i += 1
medium = substrate.copy()
grid = substrate.copy()


'''
Create a dictionary called cells that holds all the cell objects in the
environment.

Unique cell ids are assigned.
Parent ids are set to None.
x and y coordinates are randomly assigned.
The frequency of cell types is determined in accordance with the proportion
defined in the parameters section above.
'''
cells = {}
for i in range(num_cells):
    cell_id = i + 1
    x = np.random.uniform(extra['xmin'], extra['xmax'], 1)
    y = np.random.uniform(extra['ymin'], extra['ymax'], 1)

    grid_index = grid_indices[int(x),int(y)]
    grid[grid_index][cell_id] = True

    if np.random.uniform(0, 1, 1) < proportion:
        params = rjeb_params.copy()
        params['x'] = x
        params['y'] = y
        cells[cell_id] = cell.Cell(cell_id, None, params, paramup_fxn)
    else:
        params = jeb_params.copy()
        params['x'] = x
        params['y'] = y
        cells[cell_id] = cell.Cell(cell_id, None, params, paramup_fxn)

# Calculate the search radius for evaluating motility in this experiment.
search_rad = max([rjeb_params['cell_rad'],jeb_params['cell_rad']])
search_rad += max([rjeb_params['max_speed'],jeb_params['max_speed']])
search_rad = int(np.ceil(search_rad))

# Remove previous outputs
try:# csv output
    os.remove(filename)
except OSError:
    pass

fnames = os.listdir('vid')# vid folder outputs
for i in fnames:
    try:
        os.remove(str('vid/')+i)
    except OSError:
        pass

# Open a connection to the csv file that will store the summarized output.
f = open(filename, 'a')
writer = csv.writer(f, lineterminator='\n')

'''
Iterate over each cell and timepoint.
'''
current_cell_id = num_cells + 1
for timepoint in range(num_timepoints):
    # Print the current number of cells and the current timepoint.
    if (timepoint % 1) == 0:
        print('Time: '+str(timepoint+1)+' / '+str(num_timepoints))
        print('Num cells: '+str(len(cells)))

    # If the timepoint is equal to the desired time of scratching then scratch.
    if timepoint == scratch_time:
        # Identify cells within scratch bounds.
        mid_cells = grid_indices[:,scratch_left:scratch_right].flatten()
        mid_cells = [list(grid[key].keys()) for key in mid_cells]
        mid_cells = list(chain(*mid_cells))

        # Remove cells within scratch bounds.
        for i in range(len(mid_cells)):
            xpos = int(cells[mid_cells[i]].params['x'])
            ypos = int(cells[mid_cells[i]].params['y'])
            del grid[grid_indices[xpos,ypos]][mid_cells[i]]
            del cells[mid_cells[i]]

    # Determine the number of cells in the current grid column. Write to file.
    row = []
    for i in range(len(grid_indices[0])-1):
        r_out = [list(grid[key].keys()) for key in grid_indices[:,i].flatten()]
        row.append(len(list(chain(*r_out))))
    writer.writerow(row)

    # Output depiction of the epithelial sheet to a png in the vid folder.
    if (timepoint % frame_pic_mod) == 0:
        # Create the depiction. This is a density depiction (cells/element).
        img_grid = np.zeros((len(grid_indices),len(grid_indices[0])))
        for x in range(len(grid_indices)):
            for y in range(len(grid_indices[0])):
                img_grid[x][y]  += len(grid[grid_indices[x,y]])
        # Create plot.
        plt.imshow(img_grid)
        plt.savefig("vid/t%d.png" % timepoint)

    # Create a list of the cells currently in the environment.
    current_cells = set(cells.keys())
    for cell_id in list(current_cells):
        # Define the cell object.
        cell_obj = cells[cell_id]

        # Determine cell's x and y indices, grid index, substrate, and medium.
        xind, yind = int(cell_obj.params['x']), int(cell_obj.params['y'])
        grid_index = int(grid_indices[xind,yind])
        csubstrate = substrate[grid_index]
        cmedium = medium[grid_index]

        # Save the cell's old x and y coordinates.
        oldx = cell_obj.params['x']
        oldy = cell_obj.params['y']

        flags, dmatrix, dmedium = cell_obj.iterate(cell_obj.params,
                                                   csubstrate, cmedium,
                                                   extra,
                                                   grid, grid_indices,
                                                   cells,
                                                   search_rad)

        if flags['death']:
            '''
            If the cell died then remove it from grid and cells.
            '''
            del grid[grid_index][cell_id]
            del cells[cell_id]
        elif flags['division']:
            '''
            If the cell divides then create daughter cells, add them to cells
            and grid. Also, remove the parent cell from cells and grid.
            '''
            parent_id = cell_obj.cell_id
            daughter1, daughter2 = flags['division_result']

            # Create daughters.
            daughter1 = cell.Cell(current_cell_id, cell_id,
                                  daughter1,
                                  paramup_fxn)
            current_cell_id += 1
            daughter2 = cell.Cell(current_cell_id, cell_id,
                                  daughter2,
                                  paramup_fxn)
            current_cell_id += 1

            # Determine grid indices of daughters.
            d1x = int(daughter1.params['x'])
            d1y = int(daughter1.params['y'])
            d1_grid_index = grid_indices[d1x,d1y]
            d2x = int(daughter2.params['x'])
            d2y = int(daughter2.params['y'])
            d2_grid_index = grid_indices[d2x,d2y]

            # Add daughters to grid and cells.
            grid[d1_grid_index][daughter1.cell_id] = True
            grid[d2_grid_index][daughter2.cell_id] = True
            cells[daughter1.cell_id] = daughter1
            cells[daughter2.cell_id] = daughter2

            # Remove parent from grid and cells.
            del grid[grid_index][cell_id]
            del cells[cell_id]
        else:
            '''
            If the cell did not die and did not divide then handle its
            motility by removing it from grid and updating its grid position.
            '''
            del grid[grid_index][cell_id]
            xpos = int(cell_obj.params['x'])
            ypos = int(cell_obj.params['y'])
            new_grid_index = grid_indices[xpos,ypos]
            grid[new_grid_index][cell_id] = True

    num_middle_cells = grid_indices[:,scratch_left:scratch_right].flatten()
    num_middle_cells = [list(grid[key].keys()) for key in num_middle_cells]
    num_middle_cells = len(list(chain(*num_middle_cells)))
    print('Cells in middle: '+str(num_middle_cells))

# Close the connection to the csv file.
f.close()
