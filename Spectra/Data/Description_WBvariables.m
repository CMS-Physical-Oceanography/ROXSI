%% Description_WBvariables.m


        % Noah Clark
        % Created: 6/13/2023
        % CMS Summer Internship
   
        
    % Purpose: define all of the variables that are saved in WBvariables

    
    
%% China Rock Variables:

% - Bdepth: cell array containing the depth vectors recorded by China 
%            Rock's wave buoys
% - BEMEM: cell array containing the wave directions recorded by the China 
%           Rock wave buoys organized based on wave frequency and time
%          frequency and time
% - Bfreq: cell array containing the frequency range used/recorded by the
%           China Rock wave buoys
% - BGivenHsig: cell array containing the given significant wave heights
%                recorded by the China Rock wave buoys 
% - Bt_Hsig: my calculated Hsigs for each hour for the China Rock buoys
%                using the method of integrating under the frequency 
%                spectrum curve
% - Blat: a vector containing the latitudes for the 3 China Rock wave buoys
% - Blon: a vector containing the longitudes for the 3 China Rock wave buoys
% - Butm: the utm coordinate location of each of the China Rock wave buoys
% - BSee: cell array containing the wave energies recorded by the China 
%          Rock wave buoys organized based on wave frequency and time
% - Bmeanspec: the time averaged spectrum for the China Rock buoys
% - Btime: cell array containing the date-time vectors for China Rock's
%           wave buoys
% - BNormWaveDir: the incoming normal wave direction for China Rock's
%                  shoreline
% - B_TED: the total energy dissipation between each of the China Rock wave
%          buoys at each frequency (Ex: B_TED{1} --> the total energy 
%          dissipation between buoys B01 and B03) (positive value means
%          that energy has been dissipated)


%% Asilomar Variables:

% - Xdepth: cell array containing the depth vectors recorded by Asilomar's
%            wave buoys
% - XEMEM: cell array containing the wave directions recorded by the
%           Asilomar wave buoys organized based on wave frequency and time
% - Xfreq: cell array containing the frequency range used/recorded by the
%           Asilomar wave buoys
% - XGivenHsig: cell array containing the given significant wave heights
%                recorded by the China Rock wave buoys
% - Xt_Hsig: my calculated Hsigs for each hour for the Asilomar buoys
%                using the method of integrating under the frequency 
%                spectrum curve
% - Xlat: a vector containing the latitudes for the 3 Asilomar wave buoys
% - Xlon: a vector containing the longitudes for the 3 Asilomar wave buoys
% - Xutm: the utm coordinate location of each of the Asilomar wave buoys
% - XSee: cell array containing the wave energies recorded by the Asilomar 
%          wave buoys organized based on wave frequency and time
% - Xmeanspec: the time averaged spectrum for the China Rock buoys
% - Xtime: cell array containing the date-time vectors for Asilomar's
%           wave buoys
% - XavgD: the average depth measured by Asilomar's wave buoys
% - XNormWaveDir: the incoming normal wave direction for Asilomar's
%                  shoreline
% - X_TED: the total energy dissipation between each of the Asilomar wave
%          buoys at each frequency (Ex: X_TED{1} --> the total energy 
%          dissipation between buoys X01 and X03) (positive value means
%          that energy has been dissipated)


%% Model Variables:

% - M_Hsig: a vector containing the significant wave heights predicted by
%            the model 
% - ModelSee: a structure containing the "See" (wave energy) for what the 
%              model predicts for all 6 of the wave buoys; it also contains
%              the corresponding time vectors (GMT) and frequency vectors


%% Other Variables:

% - mooringLat: the latitudes of all of the instruments used in the field 
%                for the ROXSI project 
% - mooringLon: the longitudes of all of the instruments used in the field 
%                for the ROXSI project 
% - mooringtable: table of the longitudes and latitudes, planned longitudes and
%                  latitudes, easting and northing, and x and y coordinates
%                  of all of the instruments used in the field for the 
%                  ROXSI project
% - NOAA: table containing the predicted and verified tide level (NAVD;LST)
%          from NOAA Tides and Currents
% - OffWindDir: vector containing the offshore (20 km(??)) wind
%                direction
% - OffWindSpd: vector containing the offshore (20 km(??))wind speed
% - OffWindTime: vector containing the offshore (20 km(??))wind time (in
%                 Los Angeles time zone)
% - plotcolors: vector containing a list of colors to easily plot in for
%                loops
% - WindDir: a vector containing the directions of the nearshore wind
% - WindDT: a vector containing the date-times of the nearshore wind
% - WindSpd: a vector containing the speed of the nearshore wind
% - NormalDirectionDifference: a table comparing how for from the normal 
%                               direction (to the shoreline) do the incoming
%                               waves average for each wave buoy
% - EXm: a vector containing the mean bottom wave orbital excursions for
%         buoy X01 at each time (hour)
% - EXp: a vector containing the peak bottom wave orbital excursions for
%         buoy X01 at each time(hour)
% - fm: a vector containing the mean energy weighted wave frequency for
%        buoy X01 at each time (hour)
% - fp: a vector containing the peak wave frequency for buoy X01 at each
%        time (hour)
% - Tm: a vector containing the mean energy weighted wave period for
%        buoy X01 at each time (hour)
% - Tp: a vector containing the peak wave period for buoy X01 at each
%        time (hour)
% - Gfreq: a vector containing the range of frequencies used
% - Hsig: a vector containing the 6 average significant wave heights
%          corresponding to the 6 wave buoys
% - mBOV: a vector containing the mean bottom orbital velocity for buoy X01
%          at each time (hour)
% - pBOV: a vector containing the peak bottom orbital velocity for buoy X01
%          at each time (hour)
% - mCelerity: a vector containing the mean celerity of the waves for buoy
%               X01 at each time (hour)
% - pCelerity: a vector containing the peak celerity of the waves for buoy
%               X01 at each time (hour)
% - mWavelength: a vector containing the mean wavelength of the wave for
%                 buoy X01 at each time (hour)
% - pWavelength: a vector containing the peak wavelength of the wave for
%                 buoy X01 at each time (hour)
% - ShoreP: a structure containing the points which create a parallel
%            shoreline at China Rock and Asilomar
% - Table_Hsig: a table that contains the average significant wave heights
%                of each of the 6 wave buoys
% 
%










