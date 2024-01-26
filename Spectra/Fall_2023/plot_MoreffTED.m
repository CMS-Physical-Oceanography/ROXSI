%% plot_MoreffTED.m
%
% Created By: Noah Clark
% Date: 9/8/23
%
% Purpose: Determine and plot the Nielsen and observed energy dissipation
%          from the nearshore buoys/ADCPs at Asilomar and China Rock 
%          (b/t B05-B10, B10-B13, X04-X05, & X05-X11)



%% Preliminaries:

clc;clear;
load('SM&ADCP_All.mat')

%addpath('C:\Users\nsc4750\OneDrive - UNC-Wilmington\CMS Work\Summer (0516-0728)')
load('WBvariables.mat','Bmeanspec','Xmeanspec','Butm','Xutm','BSee',...
    'XSee','BEMEM','XEMEM','Xdepth','Bdepth','Btime','Xtime')

addpath('/Users/noahclark/Library/CloudStorage/OneDrive-UNC-Wilmington/CMS Work/Summer (0516-0728)/Main Functions')

a1 = 5.5; a2 = -0.2; a3 = -6.3;
kw = 1;
ff = (0.0098:0.0098:0.343);

%% From X04 to X05

% Assigning:

See1 = nanmean(XSee{3}(:,152:end)');%Xmeanspec{3}(1:36); %spotter X04
See1 = See1(1:35);
See1 = See1';
See2 = interp1(ADCP.freq,ADCP.X05.AVGspec,ff); %ADCP X05
See2 = See2';

Direc1 = meanangle(XEMEM{3}(1:35,[1:172 174:681]),2);
Direc2 = Direc1;
utm1 = Xutm{3};
utm2 = ADCP.X05.utm;
Depth1 = Xdepth{3}(152:end); %from time 21-Jun-2022 19:00 to 20-Jul-2022 4:00
Depth2 = ADCP.X05.depth(1:681); %match same time


[TED,~,~,~] = ...
    NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);

[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
TED_fer = (TED_fer1 + TED_fer2)./2;

[~,~,~,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
[~,~,~,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
TED_fep = (TED_fep1 + TED_fep2)./2;

[~,~,~,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,~,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');
TED_fem = (TED_fem1 + TED_fem2)./2;


figure(3);clf;
plot(ff,TED,'LineWidth',1.5)
hold on
plot(ff,TED_fer,'--k','LineWidth',1.5)
plot(ff,TED_fep,'--r','LineWidth',1.5)
plot(ff,TED_fem,'--b','LineWidth',1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Energy Diss')
legend('Obs TED','Nielson (f_e_R) TED','Nielson (f_e_P) TED','Nielson (f_e_M) TED')
title('Energy Dissipation from X04 to X05')




%% From X05 to X11

See1 = interp1(ADCP.freq,ADCP.X05.AVGspec,ff); %ADCP X05
See1 = See1';
See2 = interp1(ADCP.freq,ADCP.X11.AVGspec,ff); %ADCP X11
See2 = See2';

Direc1 = meanangle(XEMEM{3}(1:35,[1:172 174:681]),2);
Direc2 = Direc1;
utm1 = ADCP.X05.utm;
utm2 = ADCP.X11.utm;
Depth1 = ADCP.X05.depth(1:681); %from time 21-Jun-2022 19:00 to 20-Jul-2022 4:00
Depth2 = ADCP.X11.depth(1:681); %match same time


[TED,~,~,~] = ...
    NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);

[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
TED_fer = (TED_fer1 + TED_fer2)./2;

[~,~,~,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
[~,~,~,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
TED_fep = (TED_fep1 + TED_fep2)./2;

[~,~,~,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,~,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');
TED_fem = (TED_fem1 + TED_fem2)./2;


figure(4);clf;
plot(ff,TED,'LineWidth',1.5)
hold on
plot(ff,TED_fer,'--k','LineWidth',1.5)
plot(ff,TED_fep,'--r','LineWidth',1.5)
plot(ff,TED_fem,'--b','LineWidth',1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Energy Diss')
legend('Obs TED','Nielson (f_e_R) TED','Nielson (f_e_P) TED','Nielson (f_e_M) TED')
title('Energy Dissipation from X05 to X11')


%% From B05 to B10

% Assigning:

See1 = nanmean(BSee{3}(:,152:end)');%Xmeanspec{3}(1:36); %spotter X04
See1 = See1(1:35);
See1 = See1';
See2 = interp1(ADCP.freq,ADCP.B10.AVGspec,ff); %ADCP B10 (RIGHT NOW INCLUDING MISSING TIME FROM SEE1)
See2 = See2';

Direc1 = meanangle(BEMEM{3}(1:35,1:681),2);
Direc2 = Direc1;
utm1 = Butm{3};
utm2 = ADCP.B10.utm;
Depth1 = Bdepth{3}(213:end); %from time 21-Jun-2022 19:00 to 20-Jul-2022 4:00
Depth2 = ADCP.B10.depth(62:681); %match same time


[TED,~,~,~] = ...
    NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);

[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
TED_fer = (TED_fer1 + TED_fer2)./2;

[~,~,~,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
[~,~,~,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
TED_fep = (TED_fep1 + TED_fep2)./2;

[~,~,~,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,~,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');
TED_fem = (TED_fem1 + TED_fem2)./2;


figure(5);clf;
plot(ff,TED,'LineWidth',1.5)
hold on
plot(ff,TED_fer,'--k','LineWidth',1.5)
plot(ff,TED_fep,'--r','LineWidth',1.5)
plot(ff,TED_fem,'--b','LineWidth',1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Energy Diss')
legend('Obs TED','Nielson (f_e_R) TED','Nielson (f_e_P) TED','Nielson (f_e_M) TED')
title('Energy Dissipation from B05 to B10')


%% B10 to B13

See1 = interp1(ADCP.freq,ADCP.B10.AVGspec,ff); %ADCP B10
See1 = See1';
See2 = interp1(ADCP.freq,ADCP.B13.AVGspec,ff); %ADCP B13
See2 = See2';

Direc1 = meanangle(BEMEM{3}(1:35,62:681),2);
Direc2 = Direc1;
utm1 = ADCP.B10.utm;
utm2 = ADCP.B13.utm;
Depth1 = ADCP.B10.depth(62:681); %from time 21-Jun-2022 19:00 to 20-Jul-2022 4:00
Depth2 = ADCP.B13.depth(62:681); %match same time


[TED,~,~,~] = ...
    NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);

[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
TED_fer = (TED_fer1 + TED_fer2)./2;

[~,~,~,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
[~,~,~,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
TED_fep = (TED_fep1 + TED_fep2)./2;

[~,~,~,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,~,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');
TED_fem = (TED_fem1 + TED_fem2)./2;


figure(6);clf;
plot(ff,TED,'LineWidth',1.5)
hold on
plot(ff,TED_fer,'--k','LineWidth',1.5)
plot(ff,TED_fep,'--r','LineWidth',1.5)
plot(ff,TED_fem,'--b','LineWidth',1.5)
grid on
xlabel('Frequency (Hz)')
ylabel('Energy Diss')
legend('Obs TED','Nielson (f_e_R) TED','Nielson (f_e_P) TED','Nielson (f_e_M) TED')
title('Energy Dissipation from B10 to B13')















