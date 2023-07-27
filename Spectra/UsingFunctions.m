%% UsingFunctions.m
%CHANGE THE NAME OF SCRIPT
% Purpose:
%           To give a simple example of how the function_FricFac and
%            function_TEdis can be applied. Here, we are given a few simple
%            time averaged data sets which are all required for the
%            functions.
%
% Goal:
%           To determine and compare the observed energy dissipation(from
%            function_TEdis) and the theoretical energy dissipation (from 
%            function_FricFac).
%
%


%% Load in Given Variables

clc;clear;
load('UsingFunctions.mat')


%% Variables

%                       **Given at Start** 
%
% - See1: array containing the wave energies organized by
%          time and frequency for Point #1
%
% - See2: array containing the wave energies organized by
%          time and frequency for Point #2
%
% - Direc1: array containing the wave directions organized by
%            time and frequency for Point #1
%
% - Direc2: array containing the wave directions organized by
%            time and frequency for Point #2
%
% - ff: a vector of frequencies which See is organized by
%
% - utm1: the utm coordinates ([x y]) of the location of Point #1
%
% - utm2: the utm coordinates ([x y]) of the location of Point #2
%
% - Depth1: water depth vector as a function of time for Point #1
%
% - Depth2: water depth vector as a function of time for Point #2
% 
%
%
%                    **Determined Within Script** 
%
% - TED: the observed energy dissipation calculated using the
%         function_TEdis function
%
% - TED_fem: the theoretical energy dissipation calculated using the
%             function_FricFac function with the "m" (energy weighted
%             mean) method
%
%
%

%% Constants

%                    **Neilson's Constants** 
%
% - These are the exact constants used in the paper

              a1 = 5.5;   a2 = -0.2;    a3 = -6.3;



%                       **Roughness Length** 
%
% - In the paper, kw=0.16 but this is very dependant on the bathymetry of
%   the area
% - This given data is for an area with very rocky bathymetry, so kw is
%    estimated to be a larger value of 1

                             kw = 1;


%% Compare 

%                       **Call Functions** 
% call function_FricFac twice because it is calculated for each of the
% buoys and averaged afterwards to find the best estimate of the energy
% dissipated between the two points
TED = function_TEdis(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
[fw_m1,u_br1,fe_jm1,TED_fem1] = ...
    function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'m');
[fw_m2,u_br2,fe_jm2,TED_fem2] = ...
    function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'m');

% Here, we find the average so we can best estimate the energy dissipation
%  between the two instead of just at one or the other
TED_fem = (TED_fem1 + TED_fem2)./2;


%                         **Plotting**
figure(1);clf; set(gcf,'position',[525,470,525,525])
subplot(2,1,1)
plot(ff,See1,'-b','LineWidth',1.5)
hold on; grid on;
plot(ff,See2,'-r','LineWidth',1.5)
title('Frequency Spectrum at Point 1 and Point 2')
xlabel('Frequency (Hz)'); ylabel('Energy (m^2/Hz)');
legend(sprintf('Point 1 (Avg Depth (m) = %g\n',Depth1),...
    sprintf('Point 2 (Avg Depth (m) = %g\n',Depth2),'Location','northeast')
xlim([0 .50]);ylim([-0.2 1]);

subplot(2,1,2)
plot(ff,TED,'-m','LineWidth',2);
hold on; grid on;
xlim([0 .50]);ylim([-0.2 1]);
title('Energy Dissipation from Point 1 to Point 2');
xlabel('Frequency (Hz)');ylabel('Energy Dissipation (W/m)');

plot(ff,TED_fem,'--b','Linewidth',1.5)
legend('Obs. Dissipation','fe_M Dissipation','location','northeast')
annotation('textbox',[0.732,0.3,0.15,0.05],'String',...
    sprintf('k_w = %g\n',kw))





