%% plot_ffTED.m
% Created By: Noah Clark
% Date: 7/18/2023
%
% Purpose: To call the functions function_TEdis and function_FricFrac to 
%           determine the time-averaged observed energy dissipation and the
%           theoretical energy dissipation calculated using the friction
%           factors(representative, peak, and energy-weighted mean). The
%           observed and theoretical energy dissipation are plotted on
%           each other. This is done between each of the buoys at Asilomar
%           and China Rock.
%
%


%% Still To Do

    % Sanity check:
%   - AFTER I FINSIH, ENTER FULL ARRAYS INTO FUNCTIONS AND FIND AVERAGE
%   AFTER TO MAKE SURE THEY ARE SIMILAR


%% Preliminaries

clc;clear;
load('WBvariables.mat');

ff = (1:51).*0.0098;

a1 = 5.5; a2 = -0.2; a3 = -6.3;
kw = 1;


%% From X01 to X03:

% Assigning:
See1 = nanmean(XSee{1}(1:51,[1:172 174:end]),2);
See2 = nanmean(XSee{2}(1:51,[1:172 174:end]),2);
Direc1 = meanangle(XEMEM{1}(1:51,[1:172 174:end]),2);
Direc2 = meanangle(XEMEM{2}(1:51,[1:172 174:end]),2);
utm1 = Xutm{1};
utm2 = Xutm{2};
Depth1 = mean(Xdepth{1});
Depth2 = mean(Xdepth{2});

% Call Functions:
TED = function_TEdis(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
[~,~,fe_jr1,TED_fer1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,fe_jp1,TED_fep1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'p');
[~,~,fe_jm1,TED_fem1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,fe_jr2,TED_fer2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'r');
[~,~,fe_jp2,TED_fep2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'p');
[~,~,fe_jm2,TED_fem2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'m');

TED_fer = (TED_fer1 + TED_fer2)./2;
TED_fep = (TED_fep1 + TED_fep2)./2;
TED_fem = (TED_fem1 + TED_fem2)./2;

% Plotting:
figure(1);clf; set(gcf,'position',[0,470,525,525])
plot(ff,TED,'-m','LineWidth',2);
hold on; grid on;
xlim([0 .50]);ylim([-0.5 1]);
title('Time Averaged Energy Dissipation from X01 to X03');
xlabel('Frequency(Hz)');ylabel('Energy Dissipation [W/m]');

plot(ff,TED_fer,'--k','Linewidth',1.5)
plot(ff,TED_fep,'--r','Linewidth',1.5)
plot(ff,TED_fem,'--b','Linewidth',1.5)
yline(0)
legend('Obs. TED','fe_R TED','fe_P TED','fe_M TED','location','northeast')


%% From X03 to X04:

% Assigning:
See1 = nanmean(XSee{2}(1:51,[1:172 174:end]),2);
See2 = nanmean(XSee{3}(1:51,[1:172 174:end]),2);
Direc1 = meanangle(XEMEM{2}(1:51,[1:172 174:end]),2);
Direc2 = meanangle(XEMEM{3}(1:51,[1:172 174:end]),2);
utm1 = Xutm{2};
utm2 = Xutm{3};
Depth1 = mean(Xdepth{2});
Depth2 = mean(Xdepth{3});

% Call Functions:
TED = function_TEdis(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
[~,~,fe_jr1,TED_fer1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,fe_jp1,TED_fep1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'p');
[~,~,fe_jm1,TED_fem1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,fe_jr2,TED_fer2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'r');
[~,~,fe_jp2,TED_fep2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'p');
[~,~,fe_jm2,TED_fem2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'m');

TED_fer = (TED_fer1 + TED_fer2)./2;
TED_fep = (TED_fep1 + TED_fep2)./2;
TED_fem = (TED_fem1 + TED_fem2)./2;

% Plotting:
figure(2);clf; set(gcf,'position',[525,470,525,525])
plot(ff,TED,'-m','LineWidth',2);
hold on; grid on;
xlim([0 .50]);ylim([-0.5 1]);
title('Time Averaged Energy Dissipation from X03 to X04');
xlabel('Frequency(Hz)');ylabel('Energy Dissipation [W/m]');

plot(ff,TED_fer,'--k','Linewidth',1.5)
plot(ff,TED_fep,'--r','Linewidth',1.5)
plot(ff,TED_fem,'--b','Linewidth',1.5)
yline(0)
legend('Obs. TED','fe_R TED','fe_P TED','fe_M TED','location','northeast')


%% From B01 to B03
% HASN'T BEEN TESTED YET
%REMEMBER THAT B01 MIGHT BE DIFFERENT BECASUE ONLY 725
% Assigning:
See1 = nanmean(BSee{1}(1:51,[1:397 505:end]),2);
See2 = nanmean(BSee{2}(1:51,[1:397 505:end]),2);  % cut out the times when B01 isn't recording
Direc1 = meanangle(BEMEM{1}(1:51,[1:397 505:end]),2);
Direc2 = meanangle(BEMEM{2}(1:51,[1:397 505:end]),2);
utm1 = Butm{1};
utm2 = Butm{2};
Depth1 = mean(Bdepth{1});
Depth2 = mean(Bdepth{2});

% Call Functions
TED = function_TEdis(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
[~,~,fe_jr1,TED_fer1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,fe_jp1,TED_fep1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'p');
[~,~,fe_jm1,TED_fem1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,fe_jr2,TED_fer2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'r');
[~,~,fe_jp2,TED_fep2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'p');
[~,~,fe_jm2,TED_fem2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'m');

TED_fer = (TED_fer1 + TED_fer2)./2;
TED_fep = (TED_fep1 + TED_fep2)./2;
TED_fem = (TED_fem1 + TED_fem2)./2;

% Plotting:
figure(3);clf; set(gcf,'position',[0,70,525,525])
plot(ff,TED,'-m','LineWidth',2);
hold on; grid on;
xlim([0 .50]);ylim([-0.5 1]);
title('Time Averaged Energy Dissipation from B01 to B03');
xlabel('Frequency(Hz)');ylabel('Energy Dissipation [W/m]');

plot(ff,TED_fer,'--k','Linewidth',1.5)
plot(ff,TED_fep,'--r','Linewidth',1.5)
plot(ff,TED_fem,'--b','Linewidth',1.5)
yline(0)
legend('Obs. TED','fe_R TED','fe_P TED','fe_M TED','','location','northeast')


%% From B03 to B05

% Assigning:
See1 = nanmean(BSee{2}(1:51,:),2);
See2 = nanmean(BSee{3}(1:51,:),2);  
Direc1 = meanangle(BEMEM{2}(1:51,:),2);
Direc2 = meanangle(BEMEM{3}(1:51,:),2);
utm1 = Butm{2};
utm2 = Butm{3};
Depth1 = mean(Bdepth{2});
Depth2 = mean(Bdepth{3});

% Call Functions
TED = function_TEdis(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
[~,~,fe_jr1,TED_fer1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,fe_jp1,TED_fep1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'p');
[~,~,fe_jm1,TED_fem1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,fe_jr2,TED_fer2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'r');
[~,~,fe_jp2,TED_fep2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'p');
[~,~,fe_jm2,TED_fem2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'m');

TED_fer = (TED_fer1 + TED_fer2)./2;
TED_fep = (TED_fep1 + TED_fep2)./2;
TED_fem = (TED_fem1 + TED_fem2)./2;

% Plotting:
figure(4);clf; set(gcf,'position',[525,70,525,525])
plot(ff,TED,'-m','LineWidth',2);
hold on; grid on;
xlim([0 .50]);ylim([-0.5 1]);
title('Time Averaged Energy Dissipation from B03 to B05');
xlabel('Frequency(Hz)');ylabel('Energy Dissipation [W/m]');

plot(ff,TED_fer,'--k','Linewidth',1.5)
plot(ff,TED_fep,'--r','Linewidth',1.5)
plot(ff,TED_fem,'--b','Linewidth',1.5)
yline(0)
legend('Obs. TED','fe_R TED','fe_P TED','fe_M TED','','location','northeast')


%% From B01 to B05

% Assigning:
See1 = nanmean(BSee{1}(1:51,:),2);
See2 = nanmean(BSee{3}(1:51,:),2);  
Direc1 = meanangle(BEMEM{1}(1:51,:),2);
Direc2 = meanangle(BEMEM{3}(1:51,:),2);
utm1 = Butm{1};
utm2 = Butm{3};
Depth1 = mean(Bdepth{1});
Depth2 = mean(Bdepth{3});

% Call Functions
TED = function_TEdis(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
[~,~,fe_jr1,TED_fer1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,fe_jp1,TED_fep1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'p');
[~,~,fe_jm1,TED_fem1] = function_FricFac(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,fe_jr2,TED_fer2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'r');
[~,~,fe_jp2,TED_fep2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'p');
[~,~,fe_jm2,TED_fem2] = function_FricFac(See2,ff,Depth2,kw,a1,a2,a3,'m');

TED_fer = (TED_fer1 + TED_fer2)./2;
TED_fep = (TED_fep1 + TED_fep2)./2;
TED_fem = (TED_fem1 + TED_fem2)./2;

% Plotting:
figure(5);clf; set(gcf,'position',[900,70,525,525])
plot(ff,TED,'-m','LineWidth',2);
hold on; grid on;
xlim([0 .50]);ylim([-0.5 1]);
title('Time Averaged Energy Dissipation from B01 to B05');
xlabel('Frequency(Hz)');ylabel('Energy Dissipation [W/m]');

plot(ff,TED_fer,'--k','Linewidth',1.5)
plot(ff,TED_fep,'--r','Linewidth',1.5)
plot(ff,TED_fem,'--b','Linewidth',1.5)
yline(0)
legend('Obs. TED','fe_R TED','fe_P TED','fe_M TED','','location','northeast')


%% From MODEL B01 to B05
% Average where spectra are cocurrent in time (because this data not 
%  missing time for B01)([1:355 463:end])

%Given:
See1 = nanmean(B01_MOD.Snn(:,[1:355 463:end]),2);
See2 = nanmean(B05_MOD.Snn(:,[1:355 463:end]),2);
ffM = B01_MOD.f;
utm1 = Butm{1};
utm2 = Butm{3};
Depth1 = mean(Bdepth{1});
Depth2 = mean(Bdepth{3});

%Determining Direction:
% Not given the model directions, so I am assuming that all of the wave
% directions are perfectly straight from one buoy to the other
% Finding that direction:
Vx = utm2(1) - utm1(1);
Vy = utm2(2) - utm1(2);
V = [Vx Vy];
Length_V = sqrt(V(1)^2 + V(2)^2); %also the distance 
NorthVec = [0 1];
Length_NorthVec = 1;
dir_V = acosd(dot(NorthVec,V)/(Length_NorthVec*Length_V)) + 180;
Direc1 = dir_V.*ones(91,1);
Direc2 = Direc1;

%Calling Functions:
M_TED = function_TEdis(See1,See2,Direc1,Direc2,ffM,utm1,utm2,Depth1,Depth2);
[~,~,fe_jm1,TED_fem1] = function_FricFac(See1,ffM,Depth1,kw,a1,a2,a3,'m');
[~,~,fe_jm2,TED_fem2] = function_FricFac(See2,ffM,Depth2,kw,a1,a2,a3,'m');

TED_fem = (TED_fem1 + TED_fem2)./2;


% Plotting:
figure(6);clf; set(gcf,'position',[1390,70,525,525])
plot(ff,TED,'-m','LineWidth',2)
hold on; grid on;
plot(ffM,M_TED,'-k','LineWidth',2);
xlim([0 .3]);ylim([-0.4 0.4]);
title('Time Averaged Energy Dissipation from B01 to B05');
xlabel('Frequency(Hz)');ylabel('Energy Dissipation [W/m]');

plot(ffM,TED_fem,'--b','Linewidth',1.5)
yline(0)
legend('Obs. TED','Model TED','Neilson (f_e_M) TED','','location','northeast')


%% Saving Figures

cd 'C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)\Start Figures\TED and FF Plots'

for i = 1:6
    saveas(figure(i),sprintf('TEDs_Fig%1d.jpeg',i))
end



    %Reassign the cd to "Start 5-16"
%For the CMS Computer:
    cd 'C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)'

    
    
    
    

