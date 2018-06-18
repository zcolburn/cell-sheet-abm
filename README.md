# ABM of collective cell migration
This project is an agent-based model of collective cell migration. It consists of three files: **cell.py**, **paramup_fxn.py**, and **Main.py**. It also contains a single folder. These files and the folder are summarized below.


**cell.py** defines the *Cell* class. It holds four key pieces of data: cell id, parent cell id, cell parameters, and the parameter updater function.


**paramup_fxn.py** defines the function *paramup_fxn*, which updates cell parameters given the cell's current parameters, local substrate, local medium, and neighboring cells.


**Main.py** defines environmental variables (i.e. environment dimensions) as well as the initial conditions and perturbations.


**vid** is a folder that contains png depictions of the environment at intervals specified by the *frame_pic_mod* variable described below.


## Key parameters in Main.py
**filename** is the name of the output csv file. In this file each row is a separate timepoint and each column represents a single unit width of the grid. Each element of this matrix contains the number of cells found within that column of the environment at the corresponding timepoint.


**num_cells** is the number of cells to seed into the environment.


**num_timepoints** is the number of timepoints to simulate.


**frame_pic_mod** indicates the frequency at which environment depictions should be created in the *vid* folder.


**proportion** is the proportion of seeded cells that are initialized as rJEB (rescued junctional epidermolysis bullosa) cells. Cells not initialized as rJEB are initialized as JEB cells.


**scratch_time** is the timepoint at which scratching will occur. It is a good idea to let the sheet reach a steady state before scratching. The number of timepoints needed to reach a steady state should be determined experimentally.


**rjeb_params** and **jeb_params** are dictionaries that define initial parameters for the two cell types. Descriptions of these parameters are found in *Main.py* itself.


**x_grid_dim** and **y_grid_dim** define the horizontal and vertical dimensions of the environment, respectively.


**scratch_left** and **scratch_right** define the left and right boundaries of the scratch.
