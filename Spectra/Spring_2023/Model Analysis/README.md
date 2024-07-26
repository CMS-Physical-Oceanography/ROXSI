## ROXSI Wave Spectra (Spring 2024/Model Analysis)
**Contains all of the scripts and data where I compare the 4 models (SWAN050, SWAN050sm, COUP050, & SWAN075) to the observed data**
**Done in the Summer of 2024**

### Plotting
* **plot_AvgSpecModel.m** - Create time-averaged spectrum plots for all spotter buoys and ADCPs and for the corresponding model observation points
* **plot_LocationsModels.m** - Plot the locations of the spotter buoys and ADCPs at Asilomar and China Rock. Also, to plot the locations of the model's observation points.
* **plot_MvsO.m** - Integrate the energy per hour over the 4 frequency bands for the observed data and for the 4 models, and then plot. Determine the energy weighted mean directions and spreads for the models and for the observed in each of the 4 frequency bands. Then plot each of the modeled directions and spreads against the observed. Calculate average statistics for each model and also the statistics at each buoy in each model. Create tables for these for the energy, direction, and directional spread (used later in Taylor diagrams)
* **plot_MvsO_BT** - Determine the energy weighted mean TED, difference in direction between buoys, and difference in directional spread between the buoys for the models and for the observed in each of the 4 frequency bands. Then plot each of the modeled directions and spreads against the observed. Calculate average statistics for each model and also the statistics at each buoy in each model. Create tables for these for the TED, difference in direction, and difference in directional spread (used later in Taylor diagrams)
* **plot_TaylorDiags.m** - Use the statistics comparing the model data to the observed data, which was calculated in plot_MvsO.m and plot_MvsO_BT.m, to plot taylor diagrams for each frequency band and for each model
* **plot_TaylorDiags_Ind.m** - Use the statistics comparing the model data to the observed data, which was calculated in plot_MvsO.m and plot_MvsO_BT.m, to plot taylor diagrams of the statistics of each individual buoy for each model
* **plot_TaylorDiags_Ind_NN.m** - Plots the same as 'plot_TaylorDiags_Ind.m' but the statistics plotted on the taylor diagrams are NOT normalized (NN)
* **plot_TaylorDiags_Overall.m** - Does the same thing as 'plot_TaylorDiags.m' (?)

### Calculating
* **CombModelData.m** - Load in the provided data from the 4 models (SWAN050, SWAN050sm, COUP050, & SWAN075), add them to their own structures, and then save them in ModelResults.mat
* **calc_EFractionM.m** - Calculate the percentage of energy in the sea swell bands at each of the modeled observation points. These percentages are averaged to determine what percentage of the enrgy is in the sea and swell bands on average.
* **calc_InterpObsModel.m** - Interpolate the observed and model data to the new frequency vector. Also, cut the times to overlap each other
* **calc_pFreq.m** - Calculate the peak frequency for each hour at each buoy for the observed and model (SWAN050, SWAN050sm, COUP050, & SWAN075) data. Plot the model peak frequency results against the observed peak frequencies

### Data
* **obsVSmod_STATS.m** - Mat file containing the statistics calculated for the models and observations. Used to plot the Taylor Diagrams
