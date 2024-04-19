## ROXSI Wave Spectra (Spring 2024)
**Codes from Noah Clark to calculate wave characteristics based on wave buoy measurements from Spring 2024**
**January 2024 - May 2024**

### Plotting
* **plot_ErrBarAvgSpecs.m** - Determine 3 time periods of 1 day and 15 hours where there is only a peak in energy in the sea band, swell band, and where there is a peak at both bands. Plot the energy, wave direction, and directional spread during these times with a std envelope over each.
* **plot_EFraction.m** - Calculate and plot the fraction of total energy that is found at each day in the sea and swell bands. This is done only for buoy B03.
* **plot_DailyShoalB.m** - Calculate the observed and theoretical shoaling coefficeints at each hour and determine the EWM for each hour in the sea and swell bands. Then average the houyrs per day. This is done for all spotter buoys and ADCPs at China Rock.
* **plot_DailyShoalX.m** - Calculate the observed and theoretical shoaling coefficeints at each hour and determine the EWM for each hour in the sea and swell bands. Then average the houyrs per day. This is done for all spotter buoys and ADCPs at Asilomar.
* **plot_DailyNielsenTEDs.m** - Perform Nielsen TED analysis and determine the EWM per hour in the sea and swell bands of these TED values. Then plot these Nielsen TED daily averages. This is done for all spotter buoys and ADCPs at China Rock and Asilomar.

### Calculating
* **calc_ShoalCoeff.m** - Calculate the shoaling coefficients between the instruments by E1/E2 and compare it against the shoaling coefficients using sqrt(Cg2/Cg1) at three time periods. Those time periods were determined by isolating 3 24 hour periods where the energy was either mostly from the Sea, the Swell, or the energy was consisted of a mix similar to the entire time average

### Functions
* **bgpatch2.m** - Dr. Suanda’s function to plot an envelope of standard deviation around a series of means
* **meanangle.m** - Determines the average angle of the inputted angles. Downloaded from internet. Can input either a vector or a matrix
* **LDR.m** - Calculates an accurate estimate of only the wavelength and wave number (same as function_wavecalculateSI.m)

### Data
* **WBvariables.mat** - All of the calculated buoy data. All of data was calculated during the summer.
* **3Periods.mat** - Data calculated from script calc_ShoalCoeff.m. Data from the three selected days of where the energy was either mostly from the Sea, the Swell, or the day where the energy was consisted of a mix similar to the entire time average.
* **SM&ADCP_All.mat** - All of the data provided to me from the ADCPs and smart moorings.
* **TED.mat** - Data of the observed and Nielsen energy dissipated between all of the buoys/ADCPs at China Rock and Asilomar.
