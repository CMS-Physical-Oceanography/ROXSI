%% plot_AvgSpecModel.m
%
% Noah Clark        4/25/2024
%
% Purpose: Create time-averaged spectrums for all spotter buoys and ADCPs
%          and for the corresponding model observation points
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

clc;clear;


%% Data Organization

load('WBvariables.mat','Bmeanspec','Xmeanspec','Bfreq','Xfreq',...
    'BEMEMmean','XEMEMmean','Bdepth','Xdepth')

load('SM&ADCP_All.mat','ADCP')

addpath('Model Results')
load('ModelResults.mat')


%% China Rock

figure(1);clf;

% Spotters
plot(Bfreq{1},Bmeanspec{1},'-r','LineWidth',1.5)
hold on; grid on;
plot(Bfreq{2},Bmeanspec{2},'-b','LineWidth',1.5)
plot(Bfreq{3},Bmeanspec{3},'-g','LineWidth',1.5)

% ADCPs
plot(ADCP.freq,ADCP.B10.AVGspec,'-m','LineWidth',1.5)
plot(ADCP.freq,ADCP.B13.AVGspec,'-c','LineWidth',1.5)

% Model Results
plot(Model.B01.f,Model.B01.AVGspec,'--r','LineWidth',1.5)
plot(Model.B03.f,Model.B03.AVGspec,'--b','LineWidth',1.5)
plot(Model.B05.f,Model.B05.AVGspec,'--g','LineWidth',1.5)
plot(Model.B10.f,Model.B10.AVGspec,'--m','LineWidth',1.5)
plot(Model.B13.f,Model.B13.AVGspec,'--c','LineWidth',1.5)

title('China Rock - Time Averaged Spectra','fontsize',18)
xlim([0 0.25])
xlabel('f [Hz]','fontsize',16)
ylabel('Energy [m^2/Hz]','fontsize',16)
legend('B01 (h = 33.9m)','B03 (h = 21.0m)','B05 (h = 16.2m)',...
    'B10 (h = 9.5m)','B13 (h = 4.7m)','Model B01','Model B03',...
    'Model B05','Model B10','Model B13','NumColumns',2)


%% Asilomar


figure(2);clf;

% Spotters
plot(Xfreq{1},Xmeanspec{1},'-r','LineWidth',1.5)
hold on; grid on;
plot(Xfreq{2},Xmeanspec{2},'-b','LineWidth',1.5)
plot(Xfreq{3},Xmeanspec{3},'-g','LineWidth',1.5)

% ADCPs
plot(ADCP.freq,ADCP.X05.AVGspec,'-m','LineWidth',1.5)
plot(ADCP.freq,ADCP.X11.AVGspec,'-c','LineWidth',1.5)

% Model Results
plot(Model.X01.f,Model.X01.AVGspec,'--r','LineWidth',1.5)
plot(Model.X03.f,Model.X03.AVGspec,'--b','LineWidth',1.5)
plot(Model.X04.f,Model.X04.AVGspec,'--g','LineWidth',1.5)
plot(Model.X05.f,Model.X05.AVGspec,'--m','LineWidth',1.5)

title('Asilomar - Time Averaged Spectra','fontsize',18)
xlim([0 0.25])
xlabel('f [Hz]','fontsize',16)
ylabel('Energy [m^2/Hz]','fontsize',16)
legend('X01 (h = 20.7m)','X03 (h = 15.8m)','X04 (h = 11.2m)',...
    'X05 (h = 9.3m)','X11 (h = 2.72m)','Model X01','Model X03',...
    'Model X04','Model X05','Model X08','location','NorthEast',...
    'NumColumns',2)























