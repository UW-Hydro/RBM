# Tutorial for Running the VIC-RBM hydrologic and stream temperature model
These are the instructions for running the integrated modeling system comprised of the large-scale hydrologic model **VIC** ([Liang et al, 1994](../Documentation/References.md)), a routing model based on the work of [Lohmann et al (1996)](../Documentation/References.md) and the semi-Lagrangian water temperature model, RBM ([Yearsley, 2009,2012](../Documentation/References.md)). Model development and implementation for the large-scale hydrologic model **VIC** are described in detail on the website of the University of Washington Department of Civil and Environmental Engineering's [UW Hydro|Computational Hydrology group](https://uw-hydro.github.io/code/)
## 1. Overview
The **VIC-RBM** model is a coupled model that links three sub-models together: the **VIC** model; the routing model; and the **RBM** model. The **VIC** model is a hydrologic model that uses meteorological forcing data as input and simulates hydrologic variables such as runoff, evapotranspiration, and soil moisture, at each grid cell. The routing model takes the output from the **VIC** model along with flow network information as input, and output streamflow at specified locations along stream network. The **RBM** model takes both the streamflow results and the meteorological data as input and calculate stream temperature along stream network. Thus, the **VIC-RBM** model as a whole is able to simulate both streamflow and stream temperature at spatial scales determined by the basic **VIC** gridded network configuration. **VIC** gridded networks have been developed for the Continental United States (CONUS) at 1/16, 1/8, 1/4, 1/2 and 1 degree of latitude and longitude. As presently configured, the **VIC-RBM** model system simulates daily-averaged or subdaily stream temperatures in accordance with the timestep of input.
## 2. Model input
### 1)Input for the VIC hydrologic model
*  Meteorological forcing data at each grid cell
    *  Minimum requirement:
        - daily precipitation
        - daily maximum and minimum temperature
        - wind speed
    *  VIC includes a meteorological data disaggregator that calculates subdaily meteorological variables required by VIC from the minimum forcing variables. Details can be seen in the [VIC website](http://vic.readthedocs.io/en/develop/Documentation/Drivers/Classic/ForcingData/).
* Soil properties at each grid cell
* Vegetation information at each grid cell
### 2)Input for the routing model
*  Flow direction file – describes the topology of the river basin network.
*  Unit hydrograph file - contains the grid cell impulse response function. <br />

More details about input for Lohmann routing model can be seen [here](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/Routing/RoutingInput.shtml).
### 3)Input for the RBM stream temperature model
*  Mohseni parameters: Stream temperature as a function of time at headwaters is estimated using Mohseni method ([Mohseni et al., 1998](../Documentation/References.md)), a nonlinear regression of stream temperature on air temperature.
*  Stream channel geometry characteristics are required to calculate flow depth and velocity
using the method of [Leopold and Maddock (1953)](../Documentation/References.md).
## 3. Model output
Possible model output from the VIC-RBM model includes:

*  ***Grid-cell-based meteorological data*** at each grid cell at daily or subdaily time step (calculated
by the VIC model), including:
     * Precipitation
     * Air temperature
     * Wind speed
     * Atmospheric pressure and density
     * Vapor pressure (or vapor pressure deficit or relative humidity or specific humidity); Incoming shortwave (solar) radiation
     * Incoming longwave (or thermal) radiation
*  ***Grid-cell-based hydrologic data*** at each grid cell at daily or subdaily time step (calculated by the
VIC model), such as:
     * runoff
     * Snow cover
     * Soil moisture
     * Evapotranspiration
     * ... ...
*  ***Grid-cell-based energy data*** at each grid cell at daily or subdaily time step (calculated by the VIC model), such as:
     * Net downward shortwave radiation
     * Net downward longwave radiation
     * Net upward sensible heat
     * ... ...
*  ***Routed streamflow at specified stream locations***
     * (can be any locations within the stream network, and the resolution is the same as grid cell, i.e., stream locations are specified by indicating the grid cell in which it falls) calculated by the routing model, daily time step.
*  ***Stream temperature at specified stream locations***
     * calculated by the RBM model, at daily or subdaily time steps.
