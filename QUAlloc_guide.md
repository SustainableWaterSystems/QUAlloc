# QUAlloc guide 

This file containsa step-wise guide to install QUAlloc, general information about the model setup and how to run your first simulation.


## Installation steps

1. Clone or download the QUAlloc repository into the current working directory.

`git clone https://github.com/SustainableWaterSystems/QUAlloc`

2. QUAlloc works with Python. We recommend to install Miniconda, particularly for Python 3. Follow their instructions given at https://docs.conda.io/en/latest/miniconda.html.

3. Get the environment file from this repository [pcrglobwb_py3.yml](pcrglobwb_py3.yml) to create a conda environment where all modules required for running QUAlloc (e.g. PCRaster, netCDF4) are installed by using the following command (Note: make sure you are in the QUAlloc folder before running this line. You only need to create the environment once and it will take ~10 minutes):

`conda env create --name pcrglobwb_python3 -f pcrglobwb_py3.yml`

4. Activate the conda environment named "pcrglobwb_python3" in a command prompt:

`conda activate pcrglobwb_python3`


## QUAlloc configuration (.cfg) file

QUAlloc requires a configuration (.cfg) file to be run. This file contains all the necessary information related to your simulation (e.g. time period, study extent, etc.) and directories to your input data (e.g. climatological forcing, water management option, etc.). An example of configuration (.cfg) file is provided in the 'config' folder.

To increase the familiarity with the model setup of QUAlloc, we facilitate a self-contained example run for the Rhine-Meuse basin, providing all the required input data and .cfg file. The data can be downloaded to your local machine (228MB) from the Zenodo repository: https://doi.org/10.5281/zenodo.14511236.

Some adjustments must be made in the .cfg file before running QUAlloc:
- `inputpath =`  : set the directory where the input data is stored.
- `outputpath =` : set a directory where you prefer the output data to be reported and which you can access.
- `clone =`      : this file defines the spatial resolution and extent of your study area, and must be in the pcraster format. Some examples of clone files are given in this [repository](https://github.com/UU-Hydro/PCR-GLOBWB_model/blob/master/clone_landmask_maps/clone_landmask_examples.zip). 

Other adjustments are optional (e.g.):
- `startyear =` and `endyear =` : denote your simulation period.
- `precipitation_ncfile =` : denotes the location of your precipitation input data file, relative to `inputDir`.
- `monthly_tot =` : denote which variables you to write to netcdf at monthly resolution to `outputpath/netcdf/`. 

Some QUAlloc specific adjustments can also be made (e.g.):
- `water_quality_flag =` : set to TRUE for simulations considering the effect of surface water quality, FALSE for prescribing this capability of the model.

When adjusting input data files in your .cfg file, remember to **check your units!**


## Running QUAlloc

Ensure the correct conda environment in a command prompt: `conda activate pcrglobwb_python3`

Navigate in the command prompt to the directory where the example was downloaded:

`cd < directory_where_demo_folder_is_located/demo >`


You can start a QUAlloc run using the following command:

`python QUAlloc_model/qualloc_runner.py config/<cfg_configuration_file>.cfg`

The run time for this example (Rhine basin, one year at a monthly scale) will take 1 minute. 


## QUAlloc outputs

Output datasets are reported for data type, sector nam and source type as follows:

`<data_type>_<sector_name>_allocated_to_<source_type>_<time-step>.nc`

`<data_type>`
  - "withdrawal": refers to the water that is withdrawn at a water source level to satisfy the demands within an allocation zone
  - "demand": refers to the withdrawn water that is supplied to each location (cell) where there are demands to satisfy

`<sector_name>`
  - "domestic"
  - "irrigation"
  - "livestock"
  - "manufacture"
  - "thermoelectric"

`<source_type>`
  - "renewable_surfacewater": refers to water obtained from the surface water system components (e.g., direct runoff, base flow, interflow, etc.)
  - "renewable_groundwater": refers to water obtained from aquifers that are recharged by percolation from the upper soil layers
  - "nonrenewable_groundwater": refers to water obtained from aquifers not replenished on a human time scale
