function [WaveNumber]= function_KwavecalculateSI(period,depth)

%FUNCTION_WAVECALCULATESI calculates an accurate estimate of the wavelength, wave
%       number, and celerity of a function using Newton's Method.
%
%Input Arguments:
%period: the wave period in seconds
%height: the wave height in meters
%depth: the wave depth in meters
%Output Argument:
%Wavelength: (L) the wavelength of the wave in meters
%WaveNumber: (k=2pi/L) the wave number 
%Clerity: (c=L/T) the wave speed in meters per second

L0 = 1;
Lprev = L0;
Lnew = 0;
thresh = 0.01;
g = 9.81; %gravitational acceleration (9.81 m/s^2)
delta = 1;


while delta > thresh
    Lprev = Lnew;
    Lnew=(g*period^2)/(2*pi)*tanh((2*pi*depth)/Lprev);
    delta=abs(Lnew-Lprev);
end

% Wavelength = Lnew;
k=(2*pi)/Lnew;
WaveNumber = k;
% WaterDepthParameter = k*depth;
% WaveSteepness = height/Lnew;
% RelativeDepth = depth/Lnew;
% Clerity = Lnew/period;

end