## ROXSI Wave Spectra
**Codes from Noah Clark to calculate wave characteristics based on wave buoy measurements**

* <u>Start.m<u> - script with all code that I've completed through 6/19
* Description_WBvariables.m - a description of every variable saved in the file WBvariables.m (found in Data Folder)
  
### Plotting
* plot_WBspectro.m - makes spectrogram plots for all 6 wave buoys
* plot_X01animation.m - create and save the animation of the X01 spectrum; determine Hsig predicted by the model using the integration method; determine both peak and mean energy-weighted period, frequency, wavelength, celerity, and bottom 
* plot_WBtAvgSpec.m - plots the time averaged spectrum for all wave buoys and computes teh average depths of each buoys
* plot_WBloc.m - plot the general locations of China Rock and Asilomar and also plot the locations of the wave buoys
* plot_WBdirecto.m - plot the EMEM spectrums ("directograms") (wave direction organized based on wave frequency and time) for each wave buoy
* plot_WBcalcHsig.m - plot the calculated Hsig vs the given Hsig for the China Rock buoys; to calculate Hsig, use the integration under the time-average spectrum
* plot_WBX01TS.m - create a figure showing the time series plot of Hsig, energy-weighted mean and peak frequency, bottom orbital velocity, wavelength, and orbital excursion for buoy X01; calculate the mean and peak orbital excursion; create a 3 panel time-series figure of the spectrogram, wind speed, and tide data from NOAA for X01; make the same plot but with the datetime limits at June 27th to July 6th
  
### Calculations
* NormDir.m - determine the incoming normal wave directions for the shoreline at China Rock and for the shoreline at Asilomar; also create a table comparing if the wave directions get closer to normal as the waves near the shore
* pmBOV.m - determines the peak and mean time-weighted frequencies, periods, bottom orbital velocities, and wavelengths for all buoys at Asilomar and China Rock
* TotalED.m - calculate the values used to determine the total energy dissipation and then calculate/estimate the total energy dissipation due to wave breaking and bottom friction

### Functions
* function_NormDir.m - calculates what incoming wave angle will be normal to the shore in Monterey Bay, California
* function_wavecalculateSI.m - calculates the wavelength, wave number, and celerity 
* function_FricFac.m - determines the friction factor based on a specified frequency, buoy, and method; and a plot is created of the frequency factor curve as a function of frequency; if 0 is entered as the input frequency, the returned plot will not have a specific frequency and/or friction factor marked
* meanangle - calculates the mean of a set of angles (in degrees) based on polar considerations (Written by J.A. Dunne, 10/20/2005)
