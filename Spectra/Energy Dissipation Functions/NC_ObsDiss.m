function[TED,GeoAngle,EWMAngle,Theta] = NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,h1,h2)
%
%NC_OBSDISS determines the observed energy dissipation between two points 
%           (typically buoys) using by comparing the fluxes at the two
%           points
%
%
% Inputs:
%           - See1: array containing the wave energies organized by 
%                    time and frequency for Point #1
%           - See2: array containing the wave energies organized by 
%                    time and frequency for Point #2
%           - Direc1: array containing the wave directions organized by
%                      time and frequency for Point #1
%           - Direc2: array containing the wave directions organized by
%                      time and frequency for Point #2
%           - ff: a vector of frequencies which See is organized by
%           - utm1: the utm coordinates ([x y]) of the location of Point #1
%           - utm2: the utm coordinates ([x y]) of the location of Point #2
%           - h1: water depth vector as a function of time for Point #1
%           - h2: water depth vector as a function of time for Point #2
%
% Outputs: 
%           - TED: the total amount of energy dissipated at each time and 
%                   frequency between the two points
%           - GeoAngle: The angle made by connecting a line through the
%                        two points (reference North)
%           - EWMAngle: The energy weighted mean wave angle at each time
%           - Theta: the angle formed between the line connecting the two
%                     points and the wave propagation direction
%
%
%
% Created By: Noah Clark       Uploaded: 8/5/2023
% 


%% Preliminaries

See{1} = See1;   %arrays must all be the same size
See{2} = See2;  
Direc{1} = Direc1;
Direc{2} = Direc2;

utm{1} = utm1;
utm{2} = utm2;

h{1} = h1;
h{2} = h2;

df = ff(2) - ff(1);
T = 1./ff; %s
omega = 2.*pi.*ff;
rho = 1025; %kg/m^3    
g = 9.81; %m/s^2

loop = size(See{1}); 

H = {[],[]};
H{1} = 4.*sqrt(See{1}.*df);
H{2} = 4.*sqrt(See{2}.*df);

a = {[],[]};
a{1} = H{1}./2;
a{2} = H{2}./2;



%% Determine the Difference in Energy Flux Through Each Buoy

Energy = {[],[]};
Energy{1} = 0.5.*rho.*g.*a{1}.^2;
Energy{2} = 0.5.*rho.*g.*a{2}.^2;
Flux = {[],[]};

for xx = 1:2
    for j = 1:loop(1)
        for i = 1:loop(2) 
                % Solving linear dispersion relationship:
            Lprev = 1; Lnew = 0; thresh = 0.01; delta = 1;
            while delta > thresh
                Lprev = Lnew;
                Lnew = (9.81*T(j)^2)/(2*pi)*tanh((2*pi*h{xx}(i))/Lprev);
                delta = abs(Lnew - Lprev);
            end
            k = (2*pi)/Lnew;
            
            Cg = 0.5*(1 + (2*k*h{xx}(i))/(sinh(2*k*h{xx}(i))))...
                *(omega(j)/k);
            Flux{xx}(j,i) = Energy{xx}(j,i)*Cg;
        end
    end
end

FluxDiff = Flux{1} - Flux{2};


%% Determine delta_r


% Find the direction from one buoy to the next (relative to North)
Vx = utm{2}(1) - utm{1}(1);
Vy = utm{2}(2) - utm{1}(2);
V = [Vx Vy];
Length_V = sqrt(V(1)^2 + V(2)^2); %also the distance 

NorthVec = [0 1];
Length_NorthVec = 1;

dir_V = acosd(dot(NorthVec,V)/(Length_NorthVec*Length_V));
GeoAngle = dir_V;

WaveDir = zeros(1,loop(2));
Theta = zeros(1,loop(2));
for i = 1:loop(2)
        % Determine energy weighted mean wave direction
    xm1 = nansum(sind(Direc{1}(:,i)).*See{1}(:,i))/nansum(See{1}(:,i));
    ym1 = nansum(cosd(Direc{1}(:,i)).*See{1}(:,i))/nansum(See{1}(:,i));
    WaveDir1 = -1*(atan2d(ym1,xm1) - 90) + 180;

    xm2 = nansum(sind(Direc{2}(:,i)).*See{2}(:,i))/nansum(See{2}(:,i));
    ym2 = nansum(cosd(Direc{2}(:,i)).*See{2}(:,i))/nansum(See{2}(:,i));
    WaveDir2 = -1*(atan2d(ym2,xm2) - 90) + 180;

    Angles = [WaveDir1, WaveDir2];

    WaveDir(i) = meanangle(Angles);
    Theta(i) = WaveDir - dir_V;
end
EWMAngle = WaveDir;
delta_r = Length_V.*cosd(Theta);

% Solve for Total Energy Dissipation
TED = FluxDiff./delta_r;


end

