function[TED] = function_TEdis(See1,See2,EMEM1,EMEM2,ff,utm1,utm2,h1,h2)

%FUNCTION_TEdis determines the energy dissipation between two points 
%                 (typically buoys) at each time and frequency
%
% Inputs:
%           - See1: cell array containing the wave energies organized by 
%                    time and frequency for Point #1
%           - See2: cell array containing the wave energies organized by 
%                    time and frequency for Point #2
%           - EMEM1: cell array containing the wave directions organized by
%                     time and frequency for Point #1
%           - EMEM2:cell array containing the wave directions organized by
%                     time and frequency for Point #2
%           - ff: a vector of frequencies which See is organized by
%           - utm1: the utm coordinates ([x y]) of the location of Point #1
%           - utm2: the utm coordinates ([x y]) of the location of Point #2
%           - h1: water depth vector as a function of time for Point #1
%           - h2: water depth vector as a function of time for Point #2
%
% Outputs: 
%           - TED: the total amount of energy dissipated at each time and 
%                   frequency between the two points
%
%
% DO I NEED TO DO ANYTHING WITH ABSOLUTE VALUES SO THAT IT ISN'T NEGATIVE
% WHEN POINT 1 IS OPPOSITE OF MY POINT 1

%% Preliminaries (only for now)
%clc;clear;
%load('WBvariables.mat')
%% Manual Inputs
%Also must add the wave directions (EMEM)
%ff = (1:129).*0.0098;

% See{1} = XSee{1};   %arrays must all be the same size
% See{2} = XSee{2};
% EMEM{1} = XEMEM{1};
% EMEM{2} = XEMEM{2};
% 
% utm{1} = Xutm{1};
% utm{2} = Xutm{2};
% 
% h{1} = Xdepth{1};
% h{2} = Xdepth{2};


%% Preliminaries (REAL)

See{1} = See1;   %arrays must all be the same size
See{2} = See2;
EMEM{1} = EMEM1;
EMEM{2} = EMEM2;

utm{1} = utm1;
utm{2} = utm2;

h{1} = h1;
h{2} = h2;


df = ff(2) - ff(1);
T = 1./ff; %s
omega = 2.*pi.*ff;
rho = 1025; %kg/m^3     %density of sea water (assumption)
g = 9.81; %m/s^2

loop = size(See{1}); %doesn't matter which See is used here (bc they're same size)

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
% delta_r = L*cos(theta)

% Find the direction from one buoy to the next (relative to North)
Vx = utm{2}(1) - utm{1}(1);
Vy = utm{2}(2) - utm{1}(2);
V = [Vx Vy];
Length_V = sqrt(V(1)^2 + V(2)^2); %also the distance 

NorthVec = [0 1];
Length_NorthVec = 1;

dir_V = acosd(dot(NorthVec,V)/(Length_NorthVec*Length_V));

% Calculate delta_r
WaveDir = zeros(loop(1),loop(2));
Theta = zeros(loop(1),loop(2));
for j = 1:loop(1)
    for i = 1:loop(2)
        Angles = [EMEM{1}(j,i)-180,EMEM{2}(j,i)-180];
        WaveDir(j,i) = atand(sum(sind(Angles))/sum(cosd(Angles)));
        Theta(j,i) = WaveDir(j,i) - dir_V + 360;
    end
end

delta_r = Length_V.*cos(Theta);

% Solve for Total Energy Dissipation
TED = FluxDiff./delta_r;



%% Plotting (will not be in final function)

        % - Small energy @ time 263
        % - Medium Energy @ time 189
        % - Large energy @ time 73
        
Gfreq = (1:129).*0.0098;
figure(1);clf;
plot(Gfreq,TED(:,263),'g','LineWidth',2)
hold on
plot(Gfreq,TED(:,189),'b','LineWidth',2)
plot(Gfreq,TED(:,73),'r','LineWidth',2)
xlabel('Frequency (Hz)');ylabel('Energy Dissipation');
legend('Low Energy Time','Medium Energy Time',...
    'High Energy Time','location','northeast')
title('Total Energy Dissipation')
grid on
xlim([0 0.3])



end


