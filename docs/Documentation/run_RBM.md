# Steps to run RBM
Following describes the steps for simulating water temperature with the grid-based semi-Lagrangian model, RBM, assuming one has downloaded [**VIC** 4.2.d](https://github.com/UW-Hydro/VIC/releases/tag/VIC.4.2.d) and prepared all the necessary input files. <br />

## 1. Download and Uncompress the Files

Download VIC_RBM2.2 [HERE](http://www.hydro.washington.edu/Lettenmaier/Models/RBM/download.shtml) <br />

Unpack the compressed file “***VIC_RBM2.2.tar.gz***” and at the Unix prompt, type:
<center>**tar –xzvf VIC_RBM2.2.tar.gz**</center>

After unpacking the file, the source code and supporting files will be in the following directory:

<center>**VIC_RBM2.2**</center>

The file:
<center>**README** - a text file describing the contents of each sub-directory.</center>

The sub-directories within this directory are:

Sub-directories | Description
--- | ---
*../Perl_Scripts*  | contains pre- and post-processing scripts written in Perl
*../RBM* | - contains the Fortran 90 source and **Makefile** that builds the executable. **RBM** and places it in the directory, **../Test_Case**. Users may wish to change this by editing the makefile.
*../rout_DA*  | - contains the source code and **Makefile** for the modified routing model, rout_RBM.
*../Test_Case* | - contains the input files and executable for running the example problem, using input data for the Salmon River basin at a gridded resolution of 1/2° lat/long.
*../Tutorial* | - contains this tutorial.
*../UH_Test* | - contains the sample output from executing the routing model, **rout_RBM**, and the input data for the sample problem.
*../VIC_Forcing* | - contains the daily precipitation (mm), max/min air temperatures (°C), and wind speed (m/sec) that provides forcings for each **VIC** grid cell.
*../VIC_Input* | - contains the snowbands, soil, veg-param files and the world_veg_lib files
*../VIC_Output* | - contains the meteorologic output files, full_data_lat_long, and the hydrologic output - files, flux_lat_long, from the VIC simulations.

## 2. Create the RBM Executable
Navigate to the folder, **../RBM**, and type:
<center>**make**</center>

This will create the executable, RBM, and copy a version to the folder, **../Test_Case**, where the example problem is found.


## 3. Run the Model
### Step 1 Build the forcing function files (flow and meteorology)
In the example problem, I used **VIC4.2.d** with the two global parameter files in the directory, **run_VIC**. To create the necessary files (executable file to run VIC), you need to download [VIC4.2.d](https://github.com/UW-Hydro/VIC/releases/tag/VIC.4.2.d). The instructions to generate the executable file **vicNl** and run VIC can be found [here](http://vic.readthedocs.io/en/vic.4.2.d/). Copy the executable file **vicNl** to the folder **run_VIC**. The example global files to generate baseflow, runoff and heat budgets are also directly included in the folder **run_VIC**, which you can use directly. The two were implemented as follows:
<center>**./vicNl –g global_param_Salmon_0.5_flux - generates base flow and runoff**</center>

and
<center>**./vicNl –g global_param_Salmon_0.5_full_data - generates heats budget** </center>

In the example, the outputs from this process are copied to the folder, ***../VIC_Output***.
### Step 2 Generate topology file
In the folder, ***../Test_Case***, run the perl script, **build_network_beta.pl**, (copied from ***../Perl_Scripts***) using the direction file that was created as described in the VIC model development (see [VIC model](https://uw-hydro.github.io/code/)). All the grid cells surrounding the basin of interest must contain a negative one (-1) for purposes of determining the headwaters segments. The number of basins and sub-basins in a river system is limited only by the amount of computer memory, but there must be only a single outlet. For example, the Columbia River system could be modeled in its entirety. Modeling the Columbia River system and another river system, for example, the Fraser River, would require two separate simulations. Also, simulated river system cannot contain braided networks. The following example is from the Salmon River in Idaho:
<center>**perl build_network_beta.pl Salmon_flowd_0.5.dir Salmon_0.5.Topology**</center>

This script requires two (2) files in the command line, the direction file (**Salmon_flowd_0.5.dir**, in the example problem) and the output file, **Salmon_0.5.Topology** (in the example problem).
### Step 3 Prepare control file

Prepare a control file (see Figure 1 for an example) that describes:

Input | Note
--- | ---
Starting Date: |
Ending Data: |
Time Steps per Day:  |  This version of model simulates only daily averages, so this is always "1".
Input directory: |
Output directory: |
Topology File: |
Network File:   |       The Network File must have the suffix, \_Network.
Flow File:      |       This is the name of the file that will be created by the routing program
Heat File:      |       This is the name of the file that will be created by the routing program
Mohseni Parameters:  |  When the characters, " grid" (note the space) are missing, constant values for the Mohseni parameters are obtained from a file with the same format as the example file, **Salmon_Parameters**.
Heat Dump:      |       FALSE if there is no advective heat source (Power Plants, for example)
Ndelta:         |       Usually "2", but can be larger, particularly in the case of slower streams

The control file must have the suffix, **.Control**, and the colon (:) after the descriptive characters in each line is required. <br />
An example control file is also included in the folder **Test_case**.
### Step 4 Generate network file
Run the Perl script, **build_input.pl**, using the control file,(**Salmon_0.5.Control**) creating the stream temperature network file (**Salmon_0.5_Network** in this example) and the ordered input file names for the routing scheme. (**Rout.Cells.init** and **Rout.Cells**). See Appendix A for a description of elements in the network file. Here again, the project name (**Salmon_0.5**) is required on the command line. The suffix, **.Control**, is appended by the Perl script:

<center>**perl build_input.pl Salmon_0.5**</center>

### Step 5 Prepare routing input file
Build the forcing function files for flow and heat budget using **rout**.

Navigate to the folder, **rout_DA**, and run:
<center>**make**</center>

This will create the executable,
<center>**rout**</center>

For this example, copy **rout** to the folder **Test_Case**

Create control file, **salmon.inp_DA**, to run the modified Lohmann routing model using the file, **Rout.Cells.init**, the first time and the file, **Rout.Cells**, if the routing model, **rout** is run again for the same set of unit hydrographs. The input file (in this example, **salmon.inp_DA**) is similar to that described on the VIC model Web site with the exception of the addition the Leopold parameters. The routing model requires a file with unit hydrographs for each cell, or, in this case, a file, **UH_ALL**, with a single unit hydrograph that is applied to all grid cells.

### Step 6 Run rounting model
Run **rout** from the directory, Test_Case, using the file, **salmon.inp_DA**
<center>**./rout salmon.inp_DA**</center>

creating the direct access files for flow and heat budget,
(**Salmon_DA_flow** and **Salmon_DA_heat** in the example) . See Appendix A for a description of the state variables in the flow and heat budget files.

### Step 7 Run RBM
All necessary files have been created at this point (see Appendix A for a description of the output). Simply run the temperature model with command line files as follows:
<center>**./RBM_VIC Salmon_0.5 Salmon_0.5.Temp**</center>

where **Salmon_0.5** refers to **Salmon_0.5_Network** ("**\_Network**" is appended by the model software) and "**Salmon_0.5.Temp**" is the output file. "**Salmon_0.5.Spat**" is also created and is a file that cross-references element numbers in the output file "**Salmon_0.5.Temp**" to lat/long for
post-processing.
