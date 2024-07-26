## ROXSI Wave Spectra (Spring 2024/Model Analysis/Full Spectra Analysis)
**Make the taylor diagrams comparing the model results to the observed data for the full sea-swell band (0.04 - 0.23 Hz) instead of using 4 separate frequency bands**
**Done in the Summer of 2024**

### Plotting 
* **plot_TaylorDiags_All.m** - Use the statistics comparing the 4 models to the observed measurements (from 'Stats_FS.mat') to plot taylor diagrams. These taylor diagrams represent the models over the majority of the frequency range (0.04-0.23Hz) as opposed to breaking them up into the 4 frequency bands as done previously.

### Calculating
* **calc_MvsO_Stats_All.m** - calculate the means per hour over most of the spectra (0.04-0.23Hz) for each of the models and observed at each buoy. Then compare them by using the mod_error function to calculate the statistics that will be used to create taylor diagrams (done in 'plot_TaylorDiags_All.m')

### Data
* **Stats_FS.mat** - Contains the statistics calculated for the models and observations but calculated over one large portion of the frequency range (0.04-0.23Hz). Calculated in 'calc_MvsO_Stats_All.m'
