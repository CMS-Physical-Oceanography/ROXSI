## WBvariables.mat
**Every important variable calculated in Start.m**
* basically all other scripts load in this data at the beginning to use within the script
* in the future, the size of this file may become a concern (must stay under 25 mb)

## Dir_Spec_X01_SWAN.mat
* data from the Swan model for buoy X01 to see what the theoretical frequency spectrum would be using this model
* Equivalent to the model version of Bsee or Xsee 

## CO-OPS_9413450_met.csv
**Tide Data from NOAA**
* Units: meters
* Time: entire year of 2022
* Time-Zone: LST (Monterey Bay)
* Datum: NAVD

## MBARI2022.mat
**Offshore wind data**
* the main use is to have some type of wind data for the portion of time where the nearshore wind data was not being recorded

## mooringtable_largescale_2022.mat
**Locations of all of the sensors (not just buoys)**
* gives the lat and lon
* gives easting and northing
* gives the x and y for their created coordinate systems
* gives the planned lat and lon

## roxsi_ISPAR_wind.mat
**Wind data from a nearby area**
* Wind data from sonic anemometer at the ISPAR buoy in the ROXSI 2022 experiment. The anemometer is 4 m above the ocean surface. Measurements are taken every hour, in 27-minute bursts. Only the average is provided at the average time dtime. Wind direction is the direction the wind is blowing to, in degrees, clockwise from the geographic north (0 = north). Wind speed is in m/s. For location, look for the GPS coordinates from the power float on the ISPAR mooring
* Note: there is data missing for a portion of time where the wave buoys were collecting data

## roxsi_spotter_L2_B01_1158_reduced.mat
**Wave buoy data (B01 - furthest offshore @ China Rock) (PART 1)**
* Note: this is the data from before B01 was removed from the water on July 2nd, 2022 @00:00
* Each file contains variables: mooring ID, SN, site, latitude, X, Y, dtime, frequency, See, nfft, df, direction_nautical, list_methodsdirspec, EMEM, freqband, peak frequency, mean frequency, Hsig, and bottom depth

## roxsi_spotter_L2_B01_1150_reduced.mat
**Wave buoy data (B01 - furthest offshore @ China Rock) (PART 2)**
* Note: this is the data from after B01 was added back into the water on July 6th, 2022 @12:00
* Each file contains variables: mooring ID, SN, site, latitude, X, Y, dtime, frequency, See, nfft, df, direction_nautical, list_methodsdirspec, EMEM, freqband, peak frequency, mean frequency, Hsig, and bottom depth

## roxsi_spotter_L2_B03_1152_reduced.mat
**Wave buoy data (B03 - middle @ China Rock)**
* Each file contains variables: mooring ID, SN, site, latitude, X, Y, dtime, frequency, See, nfft, df, direction_nautical, list_methodsdirspec, EMEM, freqband, peak frequency, mean frequency, Hsig, and bottom depth

## roxsi_spotter_L2_B05_1153_reduced.mat
**Wave buoy data (B05 - most nearshore @ China Rock)**
* Each file contains variables: mooring ID, SN, site, latitude, X, Y, dtime, frequency, See, nfft, df, direction_nautical, list_methodsdirspec, EMEM, freqband, peak frequency, mean frequency, Hsig, and bottom depth

## roxsi_spotter_L2_X01_1151_reduced.mat
**Wave buoy data (X01 - furthest offshore @ Asilomar)**
* Each file contains variables: mooring ID, SN, site, latitude, X, Y, dtime, frequency, See, nfft, df, direction_nautical, list_methodsdirspec, EMEM, freqband, peak frequency, mean frequency, Hsig, and bottom depth

## roxsi_spotter_L2_X03_1157_reduced.mat
**Wave buoy data (X03 - middle @ Asilomar)**
* Each file contains variables: mooring ID, SN, site, latitude, X, Y, dtime, frequency, See, nfft, df, direction_nautical, list_methodsdirspec, EMEM, freqband, peak frequency, mean frequency, Hsig, and bottom depth

## roxsi_spotter_L2_X04_1155_reduced.mat
**Wave buoy data (X04 - closest nearshore @ Asilomar)**
* Each file contains variables: mooring ID, SN, site, latitude, X, Y, dtime, frequency, See, nfft, df, direction_nautical, list_methodsdirspec, EMEM, freqband, peak frequency, mean frequency, Hsig, and bottom depth

