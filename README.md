--------------------------------------------------------------------------------------
# CitcomVE Code and Data Processing in Python

### Contents of this repository:
- CitcomVE code (C) 
- input files
- processing codes (Python)
- some figure codes (Python)

### Running Citcom
See Shijie's tutorial (is included in this respository). 

### Current configuration:
See manuscript for details on viscosity structures, geometry, wavelengths, and periods.  

#### In the data:
Each period is denoted by "c" followed by a number between 0 and 32 (e.g. `c01`).  `c01` corresponds to a loading period of 320 Maxwell times, `c02` to 160 Maxwell times, etc (Table 2). 

Loading wavelengths are denoted as follows, '' = 5D, 'A' = 2.5D, 'B' = 10D, 'C' = 1.25D, and 'D' = 0.625D.  Box width always equals 1/2 of the wavelength.  Note that some of my python code uses the words wavelength and box width interchangably, and mostly refers to box width, not actually wavelength. 

Results are saved to folders named in the inputfile `datafile` field.  
 
##### Data file types: 
 
`period.stress.timestep`: contains instantaneous dissipation for a given period at each node at each timestep during the loading cycle. Columns are: 'stress', 'dissipation', 'viscosity'.

`velo`: contains nodes and coordinates. Columns are: 'Node','X','Z','deltaX','deltaZ','na','na'

`period.time`: contains total dissipation for the box at each timestep for a given period. Columns are: 'Timestep','Time', 'Total Dissipation','Total Elastic','Total Elastic+Dissipation','Total Work', 'Dissipation', 'Elastic Energy','Work'

`period.topo_s_timestep.dat`: contains topography for the surface of the box for each timesetp for a given period. Columns are: 'x','z','zprime','na','na'.  


### Data processing 
Code is mostly fairly well documented, and should be fairly self explanatory, but please email me if you have questions! 

### Figures
I included the code that makes the figures in our paper for examples. 

### General notes that might be useful

- If you need to change the number of timesteps, the delta_t field in the input files doesn't do anything, you have to instead change this in the Citcom.c file.  The current numbers of timesteps is pretty stable though, so I wouldn't change it unless you need to. 
- Most of the weird/annoying parts of my Python codes are because the number of timesteps changes for each period, and sometimes between wavelengths. If the scripts aren't finding files for a given case this is usually the problem. 

