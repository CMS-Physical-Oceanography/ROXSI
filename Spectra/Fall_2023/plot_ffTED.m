%% plot_ffTED.m
% Created By: Noah Clark
% Date: 7/18/2023
%
% Purpose: To call the functions NC_ObsDiss and NC_NeilsonDiss to 
%           determine the time-averaged observed energy dissipation and the
%           theoretical energy dissipation calculated using the friction
%           factors(representative, peak, and energy-weighted mean). The
%           observed and theoretical energy dissipation are plotted on
%           each other. This is done between each of the buoys at Asilomar
%           and China Rock.
%
%


%% Preliminaries

clc;clear;
load('WBvariables.mat','XSee','BSee','XEMEM','BEMEM','Xutm','Butm',...
    'Xdepth','Bdepth','B01_MOD','B05_MOD','Bmeanspec','Xmeanspec',...
    'Btime','Xtime')
load('SM&ADCP_All.mat')

%FunctionPath = dir('/Users/noahclark/Library/CloudStorage/OneDrive-UNC-Wilmington/CMS Work/Summer (0516-0728)/Main Functions');
%addpath(FunctionPath.folder)
addpath('../Summer (0516-0728)/Main Functions')

%ff = (1:51).*0.0098;
ff = (0.0098:0.0098:0.343);

a1 = 5.5; a2 = -0.2; a3 = -6.3;
kw = 1;
kw2 = 0.50;
kw3 = 0.15;
kwBeach = 0.02;


df = 0.0098;

%% From X01 to X03 (1)

% Assigning:
See1 = nanmean(XSee{1}(1:35,[1:172 174:end]),2);
See2 = nanmean(XSee{2}(1:35,[1:172 174:end]),2);
Direc1 = meanangle(XEMEM{1}(1:35,[1:172 174:end]),2);
Direc2 = meanangle(XEMEM{2}(1:35,[1:172 174:end]),2);
utm1 = Xutm{1};
utm2 = Xutm{2};
Depth1 = mean(Xdepth{1});
Depth2 = mean(Xdepth{2});

% Call Functions:
obsTED{1} = NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,~,beachTED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kwBeach,a1,a2,a3,'r');

% [~,~,fe_jp1,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
% [~,~,fe_jm1,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
[~,~,~,beachTED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kwBeach,a1,a2,a3,'r');
% [~,~,fe_jp2,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
% [~,~,fe_jm2,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');

TED_fer{1} = (TED_fer1 + TED_fer2)./2;
%TED_fep = (TED_fep1 + TED_fep2)./2;
%TED_fem = (TED_fem1 + TED_fem2)./2;

BeachTED_fer{1} = (beachTED_fer1 + beachTED_fer2)./2;

TED_fer_diff(:,1) = TED_fer{1} - obsTED{1};
%TED_fep_diff(:,1) = TED_fep - obsTED{1};
%TED_fem_diff(:,1) = TED_fem - obsTED{1};

% Plotting:
figure(1);clf; set(gcf,'position',[0,470,525,525])
plot(ff,obsTED{1},'-m','LineWidth',2);
hold on; grid on;
title('Time Averaged Energy Dissipation from X01 to X03');
xlabel('Frequency [Hz]');ylabel('Energy Dissipation [W/m]');

plot(ff,TED_fer{1},'--k','Linewidth',1.5)
% plot(ff,TED_fep,'--r','Linewidth',1.5)
% plot(ff,TED_fem,'--b','Linewidth',1.5)
yline(0)
plot(ff,BeachTED_fer{1},'--b','Linewidth',1.5)
legend('Obs. TED','Nielson TED','','Beach Nielson TED','location','northeast')
xlim([0 0.35])
ylim([0 0.8])


TEDnames(1,:) = '(1): X01 to X03';


%% From X03 to X04 (2)

% Assigning:
See1 = nanmean(XSee{2}(1:35,[1:172 174:end]),2);
See2 = nanmean(XSee{3}(1:35,[1:172 174:end]),2);
Direc1 = meanangle(XEMEM{2}(1:35,[1:172 174:end]),2);
Direc2 = meanangle(XEMEM{3}(1:35,[1:172 174:end]),2);
utm1 = Xutm{2};
utm2 = Xutm{3};
Depth1 = mean(Xdepth{2});
Depth2 = mean(Xdepth{3});

% Call Functions:
obsTED{2} = NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,~,beachTED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kwBeach,a1,a2,a3,'r');
% [~,~,fe_jp1,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
% [~,~,fe_jm1,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
[~,~,~,beachTED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kwBeach,a1,a2,a3,'r');
% [~,~,fe_jp2,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
% [~,~,fe_jm2,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');

TED_fer{2} = (TED_fer1 + TED_fer2)./2;
% TED_fep = (TED_fep1 + TED_fep2)./2;
% TED_fem = (TED_fem1 + TED_fem2)./2;

BeachTED_fer{2} = (beachTED_fer1 + beachTED_fer2)./2;

TED_fer_diff(:,2) = TED_fer{2} - obsTED{2};
% TED_fep_diff(:,2) = TED_fep - obsTED{2};
% TED_fem_diff(:,2) = TED_fem - obsTED{2};

% Plotting:
figure(2);clf; set(gcf,'position',[525,470,525,525])
plot(ff,obsTED{2},'-m','LineWidth',2);
hold on; grid on;
xlim([0 .35]);ylim([-0.1 1.3]);
title('Time Averaged Energy Dissipation from X03 to X04');
xlabel('Frequency [Hz]');ylabel('Energy Dissipation [W/m]');

plot(ff,TED_fer{2},'--k','Linewidth',1.5)
% plot(ff,TED_fep,'--r','Linewidth',1.5)
% plot(ff,TED_fem,'--b','Linewidth',1.5)
yline(0)
plot(ff,BeachTED_fer{2},'--b','Linewidth',1.5)
legend('Obs. TED','Nielson TED','','Beach Nielson TED','location','northeast')

TEDnames(2,:) = '(2): X03 to X04';


%% From X04 to X05 (3)

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


[obsTED{3},~,~,~] = ...
    NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);

[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
TED_fer{3} = (TED_fer1 + TED_fer2)./2;

[~,~,~,beachTED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kwBeach,a1,a2,a3,'r');
[~,~,~,beachTED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kwBeach,a1,a2,a3,'r');
BeachTED_fer{3} = (beachTED_fer1 + beachTED_fer2)./2;

% [~,~,~,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
% [~,~,~,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
% TED_fep = (TED_fep1 + TED_fep2)./2;
% 
% [~,~,~,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
% [~,~,~,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');
% TED_fem = (TED_fem1 + TED_fem2)./2;

TED_fer_diff(:,3) = TED_fer{3} - obsTED{3};
% TED_fep_diff(:,3) = TED_fep - obsTED{3};
% TED_fem_diff(:,3) = TED_fem - obsTED{3};

figure(3);clf;
plot(ff,obsTED{3},'m','LineWidth',1.5)
hold on
plot(ff,TED_fer{3},'--k','LineWidth',1.5)
% plot(ff,TED_fep,'--r','LineWidth',1.5)
% plot(ff,TED_fem,'--b','LineWidth',1.5)
grid on
xlabel('Frequency [Hz]')
ylabel('Energy Diss [W/m]')
yline(0)
plot(ff,BeachTED_fer{3},'--b','Linewidth',1.5)
legend('Obs. TED','Nielson TED','','Beach Nielson TED','location','northeast')
title('Energy Dissipation from X04 to X05')

TEDnames(3,:) = '(3): X04 to X05';


%% From X05 to X11 (4)

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


[obsTED{4},~,~,~] = ...
    NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);

[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
TED_fer{4} = (TED_fer1 + TED_fer2)./2;

[~,~,~,beachTED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kwBeach,a1,a2,a3,'r');
[~,~,~,beachTED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kwBeach,a1,a2,a3,'r');
BeachTED_fer{4} = (beachTED_fer1 + beachTED_fer2)./2;

% [~,~,~,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
% [~,~,~,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
% TED_fep = (TED_fep1 + TED_fep2)./2;
% 
% [~,~,~,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
% [~,~,~,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');
% TED_fem = (TED_fem1 + TED_fem2)./2;

TED_fer_diff(:,4) = TED_fer{4} - obsTED{4};
% TED_fep_diff(:,4) = TED_fep - obsTED{4};
% TED_fem_diff(:,4) = TED_fem - obsTED{4};

figure(4);clf;
plot(ff,obsTED{4},'m','LineWidth',1.5)
hold on
plot(ff,TED_fer{4},'--k','LineWidth',1.5)
% plot(ff,TED_fep,'--r','LineWidth',1.5)
% plot(ff,TED_fem,'--b','LineWidth',1.5)
grid on
xlabel('Frequency [Hz]')
ylabel('Energy Diss [W/m]')
yline(0)
plot(ff,BeachTED_fer{4},'--b','Linewidth',1.5)
legend('Obs. TED','Nielson TED','','Beach Nielson TED','location','northeast')
title('Energy Dissipation from X05 to X11')


TEDnames(4,:) = '(4): X05 to X11';



%% From B01 to B03 (5)
% HASN'T BEEN TESTED YET
%REMEMBER THAT B01 MIGHT BE DIFFERENT BECASUE ONLY 725
% Assigning:
See1 = nanmean(BSee{1}(1:35,[1:397 505:end]),2);
See2 = nanmean(BSee{2}(1:35,[1:397 505:end]),2);  % cut out the times when B01 isn't recording
Direc1 = meanangle(BEMEM{1}(1:35,[1:397 505:end]),2);
Direc2 = meanangle(BEMEM{2}(1:35,[1:397 505:end]),2);
utm1 = Butm{1};
utm2 = Butm{2};
Depth1 = mean(Bdepth{1});
Depth2 = mean(Bdepth{2});

% Call Functions
obsTED{5} = NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
% [~,~,fe_jp1,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
% [~,~,fe_jm1,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
% [~,~,fe_jp2,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
% [~,~,fe_jm2,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');

TED_fer{5} = (TED_fer1 + TED_fer2)./2;
% TED_fep = (TED_fep1 + TED_fep2)./2;
% TED_fem = (TED_fem1 + TED_fem2)./2;

[~,~,~,beachTED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kwBeach,a1,a2,a3,'r');
[~,~,~,beachTED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kwBeach,a1,a2,a3,'r');
[~,~,~,kw2_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw2,a1,a2,a3,'r');
[~,~,~,kw2_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw2,a1,a2,a3,'r');
[~,~,~,kw3_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw3,a1,a2,a3,'r');
[~,~,~,kw3_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw3,a1,a2,a3,'r');

kw2_fer{5} = (kw2_fer1 + kw2_fer2)./2;
kw3_fer{5} = (kw3_fer1 + kw3_fer2)./2;
BeachTED_fer{5} = (beachTED_fer1 + beachTED_fer2)./2;

TED_fer_diff(:,5) = TED_fer{5} - obsTED{5};
% TED_fep_diff(:,5) = TED_fep - obsTED{5};
% TED_fem_diff(:,5) = TED_fem - obsTED{5};

% Plotting:
figure(5);clf; set(gcf,'position',[0,70,525,525])
plot(ff,obsTED{5},'-m','LineWidth',2);
hold on; grid on;
xlim([0 .35]);ylim([-0.2 0.5]);
title('Time Averaged Energy Dissipation from B01 to B03');
xlabel('Frequency [Hz]');ylabel('Energy Dissipation [W/m]');

plot(ff,TED_fer{5},'--k','Linewidth',1.5)
% plot(ff,TED_fep,'--r','Linewidth',1.5)
% plot(ff,TED_fem,'--b','Linewidth',1.5)
yline(0)
plot(ff,BeachTED_fer{5},'--b','Linewidth',1.5)
legend('Obs. TED','Nielson TED','','Beach Nielson TED','location','northeast')

TEDnames(5,:) = '(5): B01 to B03';


%% From B03 to B05 (6)

% Assigning:
See1 = nanmean(BSee{2}(1:35,:),2);
See2 = nanmean(BSee{3}(1:35,:),2);  
Direc1 = meanangle(BEMEM{2}(1:35,:),2);
Direc2 = meanangle(BEMEM{3}(1:35,:),2);
utm1 = Butm{2};
utm2 = Butm{3};
Depth1 = mean(Bdepth{2});
Depth2 = mean(Bdepth{3});

% Call Functions
obsTED{6} = NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
% [~,~,fe_jp1,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
% [~,~,fe_jm1,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
% [~,~,fe_jp2,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
% [~,~,fe_jm2,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');

TED_fer{6} = (TED_fer1 + TED_fer2)./2;
% TED_fep = (TED_fep1 + TED_fep2)./2;
% TED_fem = (TED_fem1 + TED_fem2)./2;

[~,~,~,beachTED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kwBeach,a1,a2,a3,'r');
[~,~,~,beachTED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kwBeach,a1,a2,a3,'r');
[~,~,~,kw2_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw2,a1,a2,a3,'r');
[~,~,~,kw2_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw2,a1,a2,a3,'r');
[~,~,~,kw3_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw3,a1,a2,a3,'r');
[~,~,~,kw3_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw3,a1,a2,a3,'r');

kw2_fer{6} = (kw2_fer1 + kw2_fer2)./2;
kw3_fer{6} = (kw3_fer1 + kw3_fer2)./2;
BeachTED_fer{6} = (beachTED_fer1 + beachTED_fer2)./2;

TED_fer_diff(:,6) = TED_fer{6} - obsTED{6};
% TED_fep_diff(:,6) = TED_fep - obsTED{6};
% TED_fem_diff(:,6) = TED_fem - obsTED{6};

% Plotting:
figure(6);clf; set(gcf,'position',[525,70,525,525])
plot(ff,obsTED{6},'-m','LineWidth',2);
hold on; grid on;
xlim([0 .50]);ylim([-0.5 1]);
title('Time Averaged Energy Dissipation from B03 to B05');
xlabel('Frequency [Hz]');ylabel('Energy Dissipation [W/m]');

plot(ff,TED_fer{6},'--k','Linewidth',1.5)
% plot(ff,TED_fep,'--r','Linewidth',1.5)
% plot(ff,TED_fem,'--b','Linewidth',1.5)
yline(0)
plot(ff,BeachTED_fer{6},'--b','Linewidth',1.5)
legend('Obs. TED','Nielson TED','','Beach Nielson TED','location','northeast')
xlim([0 0.35])
ylim([-0.75 0.75])


TEDnames(6,:) = '(6): B03 to B05';


%% From B05 to B10 (7)

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


[obsTED{7},~,~,~] = ...
    NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);

[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
TED_fer{7} = (TED_fer1 + TED_fer2)./2;

[~,~,~,beachTED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kwBeach,a1,a2,a3,'r');
[~,~,~,beachTED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kwBeach,a1,a2,a3,'r');
[~,~,~,kw2_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw2,a1,a2,a3,'r');
[~,~,~,kw2_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw2,a1,a2,a3,'r');
[~,~,~,kw3_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw3,a1,a2,a3,'r');
[~,~,~,kw3_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw3,a1,a2,a3,'r');

kw2_fer{7} = (kw2_fer1 + kw2_fer2)./2;
kw3_fer{7} = (kw3_fer1 + kw3_fer2)./2;
BeachTED_fer{7} = (beachTED_fer1 + beachTED_fer2)./2;

% [~,~,~,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
% [~,~,~,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
% TED_fep = (TED_fep1 + TED_fep2)./2;
% 
% [~,~,~,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
% [~,~,~,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');
% TED_fem = (TED_fem1 + TED_fem2)./2;

TED_fer_diff(:,7) = TED_fer{7} - obsTED{7};
% TED_fep_diff(:,7) = TED_fep - obsTED{7};
% TED_fem_diff(:,7) = TED_fem - obsTED{7};

figure(7);clf;
plot(ff,obsTED{7},'m','LineWidth',1.5)
hold on
plot(ff,TED_fer{7},'--k','LineWidth',1.5)
% plot(ff,TED_fep,'--r','LineWidth',1.5)
% plot(ff,TED_fem,'--b','LineWidth',1.5)
grid on
xlabel('Frequency [Hz]')
ylabel('Energy Diss [W/m]')
yline(0)
plot(ff,BeachTED_fer{7},'--b','Linewidth',1.5)
legend('Obs. TED','Nielson TED','','Beach Nielson TED','location','northeast')
title('Energy Dissipation from B05 to B10')


TEDnames(7,:) = '(7): B05 to B10';


%% B10 to B13 (8)

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


[obsTED{8},~,~,~] = ...
    NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);

[~,~,~,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
[~,~,~,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
[~,~,~,kw2_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw2,a1,a2,a3,'r');
[~,~,~,kw2_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw2,a1,a2,a3,'r');
[~,~,~,kw3_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw3,a1,a2,a3,'r');
[~,~,~,kw3_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw3,a1,a2,a3,'r');

kw2_fer{8} = (kw2_fer1 + kw2_fer2)./2;
kw3_fer{8} = (kw3_fer1 + kw3_fer2)./2;
TED_fer{8} = (TED_fer1 + TED_fer2)./2;

[~,~,~,beachTED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kwBeach,a1,a2,a3,'r');
[~,~,~,beachTED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kwBeach,a1,a2,a3,'r');
BeachTED_fer{8} = (beachTED_fer1 + beachTED_fer2)./2;

% [~,~,~,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
% [~,~,~,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
% TED_fep = (TED_fep1 + TED_fep2)./2;
% 
% [~,~,~,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
% [~,~,~,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');
% TED_fem = (TED_fem1 + TED_fem2)./2;

TED_fer_diff(:,8) = TED_fer{8} - obsTED{8};
% TED_fep_diff(:,8) = TED_fep - obsTED{8};
% TED_fem_diff(:,8) = TED_fem - obsTED{8};

figure(8);clf;
plot(ff,obsTED{8},'m','LineWidth',1.5)
hold on
plot(ff,TED_fer{8},'--k','LineWidth',1.5)
% plot(ff,TED_fep,'--r','LineWidth',1.5)
% plot(ff,TED_fem,'--b','LineWidth',1.5)
grid on
xlabel('Frequency [Hz]')
ylabel('Energy Diss [W/m]')
yline(0)
plot(ff,BeachTED_fer{8},'--b','Linewidth',1.5)
legend('Obs. TED','Nielson TED','','Beach Nielson TED','location','northeast')
title('Energy Dissipation from B10 to B13')


TEDnames(8,:) = '(8): B10 to B13';



%% Save Variables

Avg_obsTED = obsTED;
Avg_NielsonTED = TED_fer;
Avg_BeachNielsonTED = BeachTED_fer;

save('TED.mat','ff','TEDnames','Avg_obsTED','Avg_NielsonTED','')

%% Save Figures

% cd '/Users/noahclark/Library/CloudStorage/OneDrive-UNC-Wilmington/CMS Work/Fall (0825)/Dec Pres New Pics'
% 
% for i = 1:8
%     saveas(figure(i),sprintf('TED_Fig%1d.jpeg',i))
% end
% 
% 
% 
% close all
% 
% % Change cd back:
% cd '/Users/noahclark/Library/CloudStorage/OneDrive-UNC-Wilmington/CMS Work/Fall (0825)'



%% B03 to B03 Observed TED (for DEC PRESENTATION)

F = figure(9);clf;

plot(ff,obsTED{6},'m','LineWidth',2)
hold on;grid on;
yline(0)
ax = gca(F);
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
xlabel('Frequency (Hz)','fontsize',22)
ylabel('Dissipation (W/m)','fontsize',22)
title('Observed Energy Dissipation from B03 to B05','fontsize',25)




%% PLOTTING WITH ONLY NIELSON TED FOR CHINA ROCK (DEC PRES)

close all
F = figure(10);clf;

    % B01 to B03
subplot(4,1,1)
plot(ff,obsTED{5},'m','Linewidth',2)
hold on; grid on;
plot(ff,TED_fer{5},'--k','Linewidth',2)
plot(ff,BeachTED_fer{5},'--r','Linewidth',2)
ax = gca(F);
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
yline(0)
legend('Observed','Nielsen Rocky Shore','Nielsen Sandy Shore','location','northeast','fontsize',14.5)
title('Energy Dissipation from B01 to B03','fontsize',23);
ylabel('Dissipation [W/m]','fontsize',17);
xlim([0 0.35])

    % B03 to B05
subplot(4,1,2)
plot(ff,obsTED{6},'m','Linewidth',2)
hold on; grid on;
plot(ff,TED_fer{6},'--k','Linewidth',2)
plot(ff,BeachTED_fer{6},'--r','Linewidth',2)
ax = gca(F);
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
yline(0)
legend('Observed','Nielsen Rocky Shore','Nielsen Sandy Shore','location','northeast','fontsize',14.5)
title('Energy Dissipation from B03 to B05','fontsize',23);
ylabel('Dissipation [W/m]','fontsize',17);
xlim([0 0.35])
ylim([-0.6 0.8])

    % B05 to B10
subplot(4,1,3)
plot(ff,obsTED{7},'m','Linewidth',2)
hold on; grid on;
plot(ff,TED_fer{7},'--k','Linewidth',2)
plot(ff,BeachTED_fer{7},'--r','Linewidth',2)
ax = gca(F);
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
yline(0)
legend('Observed','Nielsen Rocky Shore','Nielsen Sandy Shore','location','northeast','fontsize',14.5)
title('Energy Dissipation from B05 to B10','fontsize',23);
ylabel('Dissipation [W/m]','fontsize',17);
xlim([0 0.35])

    % B10 to B13
subplot(4,1,4)
plot(ff,obsTED{8},'m','Linewidth',2)
hold on; grid on;
plot(ff,TED_fer{8},'--k','Linewidth',2)
plot(ff,BeachTED_fer{8},'--r','Linewidth',2)
ax = gca(F);
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
yline(0)
legend('Observed','Nielsen Rocky Shore','Nielsen Sandy Shore','location','northeast','fontsize',14.5)
title('Energy Dissipation from B10 to B13','fontsize',23);
xlabel('Frequency [Hz]','fontsize',22);
ylabel('Dissipation [W/m]','fontsize',17);
xlim([0 0.35])


%% 

F = figure(11);clf;

    % B03 to B05
plot(ff,TED_fer{6},'--k','Linewidth',2)
hold on; grid on;
plot(ff,kw2_fer{6},'--g','Linewidth',2)
plot(ff,kw3_fer{6},'--c','Linewidth',2)
plot(ff,BeachTED_fer{6},'--r','Linewidth',2)
yline(0)
legend('Rocky Shore','k_w_2 Shore','k_w_3 Shore','Sandy Shore','location','northeast','fontsize',15)
title('Nielsen Energy Dissipation from B03 to B05','fontsize',25);
ax = gca(F);
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
ylabel('Dissipation [W/m]','fontsize',22);
xlim([0 0.35])
xlabel('Frequency [Hz]','fontsize',20);




%% 
% %% RMS Differences b/t Methods (R, P, M)
% 
% RMS_diff_fer = (1/numel(TED_fer_diff))*sum(TED_fer_diff.^2,'all');
% %RMS_diff_fer_var = var(TED_fer_diff.^2,1,'all') %TOO HIGH
% % 10.66
% %RMS_diff_fer_var = var(TED_fer_diff,1,'all') %HIGH
% % 0.743
% RMS_diff_fer_UNC = sqrt(sum((RMS_diff_fer - TED_fer_diff).^2,'all')/((numel(TED_fer_diff) - 1)*numel(TED_fer_diff)));
% 
% RMS_diff_fep = (1/numel(TED_fep_diff))*sum(TED_fep_diff.^2,'all');
% RMS_diff_fep_UNC = sqrt(sum((RMS_diff_fep - TED_fep_diff).^2,'all')/((numel(TED_fep_diff) - 1)*numel(TED_fep_diff)));
% 
% RMS_diff_fem = (1/numel(TED_fem_diff))*sum(TED_fem_diff.^2,'all');
% RMS_diff_fem_UNC = sqrt(sum((RMS_diff_fem - TED_fem_diff).^2,'all')/((numel(TED_fem_diff) - 1)*numel(TED_fem_diff)));
% 
% figure(9);clf;
% plot(1,RMS_diff_fer,'.r','MarkerSize',15)
% hold on
% errorbar(1,RMS_diff_fer,RMS_diff_fer_UNC)
% plot(2,RMS_diff_fep,'.b','MarkerSize',15)
% errorbar(2,RMS_diff_fep,RMS_diff_fep_UNC)
% plot(3,RMS_diff_fem,'.g','MarkerSize',15)
% errorbar(3,RMS_diff_fem,RMS_diff_fem_UNC)
% grid on
% xlim([0 4])
% set(gca,'xtick',[])
% ylim([0.6 0.9])
% ylabel('RMS Difference from Observed Dissipation')
% legend('Representative Method','','Peak Method','','Mean Method')
% 
% 
% %% RMS Differences per Avg Depth
% 
% for i = 1:8
%     RMS_D_fer(i) = (1/numel(TED_fer_diff(:,i)))*sum(TED_fer_diff(:,i).^2);
%     RMS_D_fep(i) = (1/numel(TED_fep_diff(:,i)))*sum(TED_fep_diff(:,i).^2);
%     RMS_D_fem(i) = (1/numel(TED_fem_diff(:,i)))*sum(TED_fem_diff(:,i).^2);
% end
% 
% %Asilomar:
% AvgDepth(1) = mean([mean(Xdepth{1}),mean(Xdepth{2})]);
% AvgDepth(2) = mean([mean(Xdepth{2}),mean(Xdepth{3})]);
% AvgDepth(3) = mean([mean(Xdepth{3}),nanmean(ADCP.X05.depth)]);
% AvgDepth(4) = mean([nanmean(ADCP.X05.depth),nanmean(ADCP.X11.depth)]);
% 
% %China Rock:
% AvgDepth(5) = mean([mean(Bdepth{1}),mean(Bdepth{2})]);
% AvgDepth(6) = mean([mean(Bdepth{2}),mean(Bdepth{3})]);
% AvgDepth(7) = mean([mean(Bdepth{3}),nanmean(ADCP.B10.depth)]);
% AvgDepth(8) = mean([nanmean(ADCP.B10.depth),nanmean(ADCP.B13.depth)]);
% 
% 
% figure(10);clf;
% plot(AvgDepth,RMS_D_fer,'.r','MarkerSize',13)
% hold on
% plot(AvgDepth,RMS_D_fep,'.b','MarkerSize',13)
% plot(AvgDepth,RMS_D_fem,'.g','MarkerSize',13)
% xlabel('Average Depth b/t Instruments (m)')
% ylabel('RMS Difference')
% legend('Representative Method','Peak Method','Mean Method')
% title('RMS Difference b/t Obs and Nielson TED')
% xlim([5 28])
% grid on
% 




%%
% %% From B01 to B05
% 
% % Assigning:
% See1 = nanmean(BSee{1}(1:51,:),2);
% See2 = nanmean(BSee{3}(1:51,:),2);  
% Direc1 = meanangle(BEMEM{1}(1:51,:),2);
% Direc2 = meanangle(BEMEM{3}(1:51,:),2);
% utm1 = Butm{1};
% utm2 = Butm{3};
% Depth1 = mean(Bdepth{1});
% Depth2 = mean(Bdepth{3});
% 
% % Call Functions
% TED = NC_ObsDiss(See1,See2,Direc1,Direc2,ff,utm1,utm2,Depth1,Depth2);
% [~,~,fe_jr1,TED_fer1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'r');
% [~,~,fe_jp1,TED_fep1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'p');
% [~,~,fe_jm1,TED_fem1] = NC_NeilsonDiss(See1,ff,Depth1,kw,a1,a2,a3,'m');
% [~,~,fe_jr2,TED_fer2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'r');
% [~,~,fe_jp2,TED_fep2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'p');
% [~,~,fe_jm2,TED_fem2] = NC_NeilsonDiss(See2,ff,Depth2,kw,a1,a2,a3,'m');
% 
% TED_fer = (TED_fer1 + TED_fer2)./2;
% TED_fep = (TED_fep1 + TED_fep2)./2;
% TED_fem = (TED_fem1 + TED_fem2)./2;
% 
% % Plotting:
% figure(7);clf; set(gcf,'position',[900,70,525,525])
% plot(ff,TED,'-m','LineWidth',2);
% hold on; grid on;
% xlim([0 .50]);ylim([-0.5 1]);
% title('Time Averaged Energy Dissipation from B01 to B05');
% xlabel('Frequency(Hz)');ylabel('Energy Dissipation [W/m]');
% 
% plot(ff,TED_fer,'--k','Linewidth',1.5)
% plot(ff,TED_fep,'--r','Linewidth',1.5)
% plot(ff,TED_fem,'--b','Linewidth',1.5)
% yline(0)
% legend('Obs. TED','fe_R TED','fe_P TED','fe_M TED','','location','northeast')
% 
% 
% %% From MODEL B01 to B05
% % Average where spectra are cocurrent in time (because this data not 
% %  missing time for B01)([1:355 463:end])
% 
% %Given:
% See1 = nanmean(B01_MOD.Snn(:,[1:355 463:end]),2);
% See2 = nanmean(B05_MOD.Snn(:,[1:355 463:end]),2);
% ffM = B01_MOD.f;
% utm1 = Butm{1};
% utm2 = Butm{3};
% Depth1 = mean(Bdepth{1});
% Depth2 = mean(Bdepth{3});
% 
% %Determining Direction:
% % Not given the model directions, so I am assuming that all of the wave
% % directions are perfectly straight from one buoy to the other
% % Finding that direction:
% Vx = utm2(1) - utm1(1);
% Vy = utm2(2) - utm1(2);
% V = [Vx Vy];
% Length_V = sqrt(V(1)^2 + V(2)^2); %also the distance 
% NorthVec = [0 1];
% Length_NorthVec = 1;
% dir_V = acosd(dot(NorthVec,V)/(Length_NorthVec*Length_V)) + 180;
% Direc1 = dir_V.*ones(91,1);
% Direc2 = Direc1;
% 
% %Calling Functions:
% M_TED = NC_ObsDiss(See1,See2,Direc1,Direc2,ffM,utm1,utm2,Depth1,Depth2);
% [~,~,fe_jm1,TED_fem1] = NC_NeilsonDiss(See1,ffM,Depth1,kw,a1,a2,a3,'m');
% [~,~,fe_jm2,TED_fem2] = NC_NeilsonDiss(See2,ffM,Depth2,kw,a1,a2,a3,'m');
% 
% TED_fem = (TED_fem1 + TED_fem2)./2;
% 
% 
% % Plotting:
% figure(8);clf; set(gcf,'position',[1390,70,525,525])
% plot(ff,TED,'-m','LineWidth',2)
% hold on; grid on;
% plot(ffM,M_TED,'-k','LineWidth',2);
% xlim([0 .3]);ylim([-0.4 0.4]);
% title('Time Averaged Energy Dissipation from B01 to B05');
% xlabel('Frequency(Hz)');ylabel('Energy Dissipation [W/m]');
% 
% plot(ffM,TED_fem,'--b','Linewidth',1.5)
% yline(0)
% legend('Obs. TED','Model TED','Neilson (f_e_M) TED','','location','northeast')






%% Saving Figures


%cd '/Users/noahclark/Library/CloudStorage/OneDrive-UNC-Wilmington/CMS Work/Important Images/ffTEDplots'
 
% for i = 1:8
%     saveas(figure(i),sprintf('TEDs_Fig%1d.jpeg',i))
% end
% 
% 
% 
%     %Reassign the cd to "Start 5-16"
% %For the CMS Computer:
%     cd 'C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)'

    
    
    
    

