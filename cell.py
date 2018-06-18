# -*- coding: utf-8 -*-
"""
Created on Sun May 27 19:54:25 2018

@author: ZacharyColburn
"""

class Cell(object):
    '''
    An object to hold essential cell data (location, unique id, parent id, and
    parameter updater function).
    '''
    
    def __init__(self, 
                                 
                 cell_id, # Cell id
                 parent_id, # Parent cell id
                 
                 params, # Parameters
                 paramup_fxn): # Parameter updater function
        '''
        Initialize a cell object. Specify a unique id number, its parent id 
        number (None if it has no parent), a dictionary of cell parameters, 
        and a function for updating paramters at each timepoint given matrix 
        data, medium data, and neighboring cell data.
        '''
        self.cell_id = cell_id
        self.parent_id = parent_id
        self.params = params
        self.paramup_fxn = paramup_fxn
        
        
    def iterate(self, 
                params, # Dictionary of parameters for the current cell.
                matrix, # Dictionary of matrix values, key = grid index.
                medium, # Dictionary of medium values, key = grid index.
                extra, # Dictionary of environment dimenions.
                grid, # Dictionary specifying the cells in each grid index.
                grid_indices, # Matrix for looking up grid index given x, y
                cells, # Dictionary of cell objects in the environment
                search_rad): # Search radius = max cell speed + max cell radius
        '''
        Evaluate the behavior of a single cell at a single timepoint.
        
        Updates cell parameters.
        
        Returns a dictionary of flags (flags) indicating cell death, changes 
        in the matrix composition (dmatrix), and changes in medium composition 
        (dmedium).
        '''
        flags, params, dmatrix, dmedium = self.paramup_fxn(params, 
                                                           matrix, 
                                                           medium, 
                                                           extra, 
                                                           grid, 
                                                           grid_indices, 
                                                           cells, 
                                                           search_rad)
        self.params = params # Update cell parameters
        return flags, dmatrix, dmedium
