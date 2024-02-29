## ROXSI Wave Spectra (Fall 2023)
**Codes from Noah Clark to calculate wave characteristics based on wave buoy measurements from Fall 2023**
**August 2023 - December 2023**

### Plotting
* **plot_CompObsNiel.m** - Plot the integrated obs diss vs the integrated nielsen diss. These are integrated over the sea and swell bands separately
* **plot_DirSpread.m** - Plot the time-averaged mean directional spread for all spotter buoys
* **plot_ffTED.m** - Call the functions NC_ObsDiss.m and NC_NeilsonDiss.m to determine the time-averaged observed energy dissipation and the theoretical energy dissipation. The observed and theoretical energy dissipation are plotted on top of each other. This is done between each of the buoys at Asilomar and China Rock
* **plot_meanDir.m** - Determine and plot the mean wave direction spectrums at buoy B01 and X01
* **plot_MoreAvgSpec.m** - Create time-averaged spectrums for all spotter buoys and ADCPs. Also plot the time-averaged directions in one figure
* **plot_MoreffTED.m** - Determine and plot the Nielsen and observed energy dissipation from the nearshore buoys/ADCPs at Asilomar and China Rock (b/t B05-B10, B10-B13, X04-X05, & X05-X11)
* **plot_MoreLocations.m** - Plot the locations of the smart moorings and ADCPs
* **plot_NonAvgTED.m** - Determine and plot the Nielsen and observed dissipation using the functions NC_NeilsonDiss.m and NC_ObsDiss.m at each time and frequency. This is done between all buoys and ADCPs to determine how similar the dissipations are at all points in time and frequencies by plotting spectrograms of these dissipations.
* **plot_Vector_dir.m** - Plot the positions of the China Rock buoys in utm coordinates and then plot the vectors of the true wave direction in one color and then the direction that is accounted for by Snell’s Law as a vector in another color
These Snell’s directions are calculated in the script calc_Snells.m

### Calculations
* **calc_Snells.m** - Use Snell's Law to determine the angle of refraction in the Sea and Swell bands. These angles are saved in the mat file Snells.mat

### Functions
* **LDR.m** - Calculates an accurate estimate of only the wavelength and wave number (same as function_wavecalculateSI.m)
**meanangle.m** - Determines the average angle of the inputted angles. Downloaded from internet. Can input either a vector or a matrix

### Data
* **WBvariables.mat** - All of the calculated buy data. All of the data was calculated over the summer 2023.

* **Snells.mat** - The true observed directions (average) and the calculated snells directions (average) at buoys B03 and B05. These directions are separated into the Sea and Swell bands 
* **TED.mat** - Data of the observed and Nielsen energy dissipated between all of the buoys/ADCPs at China Rock and Asilomar
