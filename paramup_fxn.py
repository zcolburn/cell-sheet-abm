# -*- coding: utf-8 -*-
"""
Created on Sun May 27 19:54:25 2018

@author: ZacharyColburn
"""

import numpy as np
from itertools import chain
from math import pi

def paramup_fxn(params, 
                csubstrate, 
                cmedium, 
                extra, 
                grid, 
                grid_indices, 
                cells, 
                search_rad):
    '''
    Function for updating cell, matrix, and medium parameters.
    
    Returns the dictionaries flags, params, dmatrix, dmedium.
    '''
    # Create empty dictionaries.
    flags = {}
    dmatrix = {}
    dmedium = {}
    
    # Get the truncated cell coordinates.
    xind, yind = int(params['x']), int(params['y'])
    
    '''
    Determine boundaries for identifying neighbors then identify neighbor ids.
    
    nlx = neighbor lower x boundary
    nux = neighbor upper x boundary
    nly = neighbor lower y boundary
    nuy = neighbor upper y boundary
    
    neighbors = list of neighbor cell ids
    '''
    nlx = int(xind-search_rad) if (xind-search_rad)  >= 0 else 0
    
    nux = (int(np.ceil(xind+search_rad)) 
        if (xind-search_rad)  < len(grid_indices) 
        else len(grid_indices)-1)
    
    nly = int(yind-search_rad) if (yind-search_rad)  >= 0 else 0
    nuy = (int(np.ceil(yind+search_rad)) 
        if (yind-search_rad)  < len(grid_indices[0]) 
        else len(grid_indices[0])-1)
    
    neighbors = grid_indices[nlx:nux,nly:nuy].flatten()
    neighbors = [list(grid[key].keys()) for key in neighbors]
    neighbors = list(chain(*neighbors))
    
    '''
    Check to see if the cell dies. If it does then cease execution and return 
    the flags, params, dmatrix, and dmedium dictionaries.
    '''
    flags['death'] = False
    if np.random.uniform(0, 1, 1) < params['death_rate']:
        flags['death'] = True
        return flags, params, dmatrix, dmedium
    
    '''
    Check to see if the cell divides. If it does then create daughter cells 
    and return the flags, params, dmatrix, and dmedium dictionaries.
    '''
    flags['division'] = False
    if np.random.uniform(0, 1, 1) <= params['division_rate']:
        flags['division'] = True
        
        # Calculate the separation dist required for division by this cell.
        min_separation_dist = params['min_division_dist']
        separation_dist = min_separation_dist + np.random.uniform(0, 1, 1)
        
        # Get the number of orientations to attempt for division.
        steps_to_test = params['div_steps']
        
        # Set px and py to the cell's x and y coordinates.
        px = params['x']
        py = params['y']
        
        '''
        Create new daughter cells and modify them using the divider function 
        below.
        '''
        daughter1 = params.copy()
        daughter2 = params.copy()
        
        def divider(d1, d2):
            return d1, d2
        daughter1, daughter2 = divider(daughter1, daughter2)
        
        '''
        Randomly attempt to divide in different orientations. Choose the 
        orientation that maximizes the distances to neighbors.
        
        The variable one_is_valid will be True if any of the division 
        iterations are possible. If none are possible then the cell won't 
        divide.
        '''
        outer_min_dist = 0
        one_is_valid = False
        for step in range(steps_to_test):
            '''
            Generate daughter coordinates by randomly choosing an x value 
            between -0.5*separation_dist and 0.5*separation_dist, then choose 
            a y value (randomly above or below the cell's y value) such that 
            the distance between the daughters is equal to the separation 
            distance defined above.
            '''
            d_x = (np.random.uniform(0, 1, 1)-0.5)*separation_dist
            y_sign = 1 if np.random.uniform(0, 1, 1) > 0.5 else -1
            d_y = np.sqrt(d_x**2+(0.5*separation_dist)**2)*y_sign
            d1_x = d_x + px
            d1_y = d_y + py
            d2_x = -d_x + px
            d2_y = -d_y + py
            
            '''
            Check whether the daughter cells are located within the 
            environment.
            '''
            if not (((d1_x >= extra['xmin']) and (d2_x >= extra['xmin'])) and 
                ((d1_x < extra['xmax']) and (d2_x < extra['xmax'])) and
                ((d1_y >= extra['ymin']) and (d2_y >= extra['ymin'])) and 
                ((d1_y < extra['ymax']) and (d2_y < extra['ymax']))):
                continue # Skip to next iteration.
            
            '''
            If these are the first set of daughter coordinates to lie within 
            the environment, then set the best daughter coordinates to the 
            current daughter coordinates.
            '''
            if not one_is_valid:
                bd1_x, bd1_y, bd2_x, bd2_y = d1_x, d1_y, d2_x, d2_y
            one_is_valid = True
            
            # If there are no neighbors then use the current coordinates.
            if len(neighbors) == 0:
                break
            if len(neighbors) > 0:
                '''
                If there are neighbors then determine their distance from the 
                daughter cells. Keep track of the shortest distance from each 
                daughter and its nearest neighbors. For the current 
                orientation store the value in step_min_dist. Update 
                outer_min_dist after checking all the neighbors.
                '''
                step_min_dist = np.Inf
                for key in neighbors:
                    val = cells[key]
                    di_1 = np.sqrt((
                            np.square(d1_x-val.params['x'])+
                            np.square(d1_y-val.params['y'])))
                    di_2 = np.sqrt((
                            np.square(d2_x-val.params['x'])+
                            np.square(d2_y-val.params['y'])))
                    step_min_dist = min([step_min_dist, di_1, di_2])
                if step_min_dist > outer_min_dist:
                    bd1_x, bd1_y, bd2_x, bd2_y = d1_x, d1_y, d2_x, d2_y
                    outer_min_dist = step_min_dist
        
        '''
        If the cell is going to divide then return the appropriate values, 
        otherwise continue evaluating.
        '''
        if (outer_min_dist > min_separation_dist) and one_is_valid:
            daughter1['x'] = bd1_x
            daughter1['y'] = bd1_y
            daughter2['x'] = bd2_x
            daughter2['y'] = bd2_y
            flags['division_result'] = (daughter1, daughter2)
            return flags, params, dmatrix, dmedium
        else:
            flags['division'] = False
    
    '''
    Calculate matrix adhesion.
    '''
    mat_ad = 1 # np.mean([params['b4_integrin'], params['paxillin']])
    
    '''
    Calculate cell-cell adhesion.
    '''
    
    
    '''
    Calculate cell polarity.
    '''
    params['prev_polarity'] = (params['polarity'] 
    if 'polarity' in params.keys() 
    else np.random.uniform(0,1,1))
    
    randnum = np.random.uniform(0,1,1)
    params['polarity'] = (randnum+params['b4_integrin']) / 2
    
    '''
    Determine cell motility.
    
    If this is the first timepoint then randomly set xvi and yvi. xvi and yvi 
    are the displacement vectors for cell motility.
    '''
    if ('xvi' not in params.keys()) and ('yvi' not in params.keys()):
        nums = np.random.uniform(0, params['max_speed'], 2)
        params['xvi'], params['yvi'] = nums
    
    '''
    Update xvi and yvi.
    '''
    xvi = ((params['xvi']*params['polarity']+
           np.random.uniform(-params['mov_chg'],params['mov_chg'],1))*mat_ad)
    yvi = ((params['yvi']*params['polarity']+
           np.random.uniform(-params['mov_chg'],params['mov_chg'],1))*mat_ad)
           
    # If the cell speed is greater than the max speed then shorten the vector.
    speed = np.sqrt(np.square(xvi)+np.square(yvi))
    if speed >= params['max_speed']:
        speed_correcter = params['max_speed']/speed
        xvi *= speed_correcter
        yvi *= speed_correcter
    
    '''
    Define the coordinates for the cell (cx, cy) and destination (mx, my).
    '''
    cx, cy = params['x'], params['y']
    pc = np.array([cx, cy])
    mx, my = np.asscalar(cx + xvi), np.asscalar(cy + yvi)
    pm = np.array([mx, my])
    def get_term_pos(cx,cy,mx,my,nx,ny,ld,crad,nrad,cmd):
        '''
        A function to calculate the cell's destination given possible 
        collisions with neighbors.
        
        cx = cell's x
        cy = cell's y
        mx = destination x
        my = destination y
        nx = neighbor x
        ny = neighbor y
        ld = distance of neighbor to the line between the cell and destination
        crad = cell radius
        nrad = neighbor radius
        cmd = cell movement distance
        '''
        # Ensure that all x and y coordinates are scalars.
        cx = np.asscalar(np.array(cx))
        cy = np.asscalar(np.array(cy))
        nx = np.asscalar(np.array(nx))
        ny = np.asscalar(np.array(ny))
        
        '''
        If the sum of cell and neighbor radii is less than the distance 
        between the neighbor and the line connecting the cell and destination, 
        then the path to the destination is clear (at least with respect to 
        this neighbor). Consequently, destination x, destination y, and cmd
        should be returned.
        '''
        trad = crad + nrad
        if ld >= trad:
            return mx, my, cmd
        else:
            '''
            Determine the angle between neighbor, cell, and destination where
            cell is the vertex.
            
            If the angle is larger than pi/3 then allow cell to move to its 
            destination.
            '''
            a = np.array([mx,my])
            b = np.array([cx,cy])
            c = np.array([nx,ny])
            ba = a - b
            bc = c - b
            
            denominator = np.linalg.norm(ba)*np.linalg.norm(bc)
            if np.abs(denominator) < 0.0001:
                cos_angle = 0
            else:
                cos_angle = np.dot(ba, bc)/denominator
            angle = np.arccos(cos_angle)
            
            if angle > (pi/3):
                return mx, my, cmd
            
            '''
            Define the line connecting the cell and destination.
            
            If the slope of the line is infinite then set it to 1000.
            '''
            denominator = cx-nx
            if denominator == 0:
                m = 10**3
            else:
                m = (cy-ny)/denominator
            b = cy-m*cx
            
            '''
            Determine the collision coordinates.
            
            There can be 0, 1, or 2 collisions.
            
            If there are no collisions then return the destination coordinates.
            '''
            a = nx
            d = ny
            e = crad
            f = nrad
            
            coeff_a = 1+m**2
            coeff_b = -2*a-2*d*m+2*b*m
            coeff_c = a**2+d**2-2*d*b+b**2-e-f
            quadratic_params = ([np.asscalar(np.array(coeff_a)), 
                                 np.asscalar(np.array(coeff_b)), 
                                 np.asscalar(np.array(coeff_c))])
            roots = np.roots(quadratic_params)
            
            roots = roots[np.isreal(roots)]
            if len(roots) == 0:# No collisions
                return mx, my, cmd
            elif len(roots) == 1:# 1 collision
                # Set x and y for the collision.
                x = roots[0]
                y = m * x + b
                '''
                If the collision occurs between the cell and destination then 
                return the new x, new y, and new cmd. Otherwise, return the 
                destination and original cmd.
                '''
                if ((((x >= cx) and (x <= mx)) or 
                     ((x <= cx) and (x >= mx))) and 
                    (((y >= cy) and (y <= my)) or 
                     ((y <= cy) and (y >= my)))):
                    cmd = np.sqrt(np.square(cx-x)+np.square(cy-y))
                    return x, y, cmd
                else:
                    return mx, my, cmd
            else:# 2 collisions
                '''
                Determine the x and y coordinates for the two collisions. 
                Calculate the cmd for each possible collision.
                '''
                x1 = roots[0]
                y1 = m * x1 + b
                x2 = roots[1]
                y2 = m * x2 + b
                cmd1 = np.sqrt(np.square(cx-x1)+np.square(cy-y1))
                cmd2 = np.sqrt(np.square(cx-x2)+np.square(cy-y2))
                
                '''
                Check whether each collision occurs between the cell and 
                destination.
                
                spot1 = True/False whether collision is between these points.
                spot2 = True/False whether collision is between these points.
                '''
                spot1 = False
                spot2 = False
                if ((((x1 >= cx) and (x1 <= mx)) or 
                     ((x1 <= cx) and (x1 >= mx))) and 
                    (((y1 >= cy) and (y1 <= my)) or 
                     ((y1 <= cy) and (y1 >= my)))):
                    spot1 = True
                elif ((((x2 >= cx) and (x2 <= mx)) or 
                       ((x2 <= cx) and (x2 >= mx))) and 
                    (((y2 >= cy) and (y2 <= my)) or 
                     ((y2 <= cy) and (y2 >= my)))):
                    spot2 = True
                    
                '''
                Return the coordinates of the first collision.
                '''
                if ((spot1 and not spot2) or 
                    (spot1 and spot2 and (cmd1 < cmd2))):
                    return x1, y1, cmd1
                elif ((spot2 and not spot1) or 
                    (spot2 and spot1 and (cmd2 < cmd1))):
                    return x2, y2, cmd2
            # Just in case, return the destination nothing has been returned.
            return mx, my, cmd
    
    # Set movement to the default destination.
    new_x, new_y, cmd_min = mx, my, np.sqrt(np.square(cx-mx)+np.square(cy-my))
    if len(neighbors) > 0:
        '''
        If there are neighbors then evaluate possible collisions.
        '''
        cmd = np.sqrt(np.square(cx-mx)+np.square(cy-my))# Maximum possible cmd
        for key in neighbors:
            val = cells[key]# val is the neighbor cell object.
            nx, ny = val.params['x'], val.params['y']# Neighbor coordinates
            pn = np.array([nx, ny])# Neighbor coordinates array
            # Get distance of neighbor to line connecting cell and destination.
            p2mp1 = [np.asscalar(pm[0]-pc[0]),np.asscalar(pm[1]-pc[1])]
            p3mp1 = [np.asscalar(pc[0]-pn[0]),np.asscalar(pc[1]-pn[1])]
            ld = np.abs(np.cross(p2mp1, p3mp1))/np.linalg.norm(p2mp1)
            
            # Get collision x, y, and cmd values.
            t_x, t_y, t_cmd = get_term_pos(cx,cy,
                                           mx,my,
                                           nx,ny,
                                           ld,
                                           params['cell_rad'],
                                           val.params['cell_rad'],
                                           cmd)
            
            '''
            If the collision distance is less than the previously determined 
            distance then update the cell's new coordinates and cmd.
            '''
            if t_cmd < cmd_min:
                new_x, new_y, cmd_min = t_x, t_y, t_cmd
    
    '''
    If the destination is within the bounds of the environment, then update 
    the x and y coordinates. Update the x and y vectors. If the destination is 
    outside the bounds of the environment then 1) do not update the x and y 
    coordinates and 2) set the x and y vectors to 0.
    '''
    if (((new_x >= extra['xmin']) and (new_y >= extra['ymin'])) and 
                ((new_x < extra['xmax']) and (new_y < extra['ymax']))):
        params['xvi'] = new_x - params['x']
        params['yvi'] = new_y - params['y']
        params['x'] = new_x
        params['y'] = new_y
    else:
        params['xvi'] = 0
        params['yvi'] = 0
        
    
    # Return the flags, params, dmatrix, and dmedium dictionaries.
    return flags, params, dmatrix, dmedium
    