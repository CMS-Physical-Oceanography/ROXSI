%% plot_MoreLocations.m
%
% Noah Clark
% Created: 9/1/2023
%
% Purpose: To plot the locations of the smart moorings and the ADCPs


clc;clear;

%% Load in Data

%addpath('C:\Users\nsc4750\OneDrive - UNC-Wilmington\CMS Work\Fall (0825)\ADCP Data')
addpath('ADCP Data')
ADCPnames = {'roxsi_signature_L2_B10_See.mat',...
    'roxsi_signature_L2_B13_See.mat','roxsi_signature_L2_X05_See.mat',...
    'roxsi_signature_L2_X11_See.mat'};

% for i = 1:length(ADCPnames)
%     load(ADCPnames{i});
%     Alat{i} = A.latitude;
%     Alon{i} = A.longitude;
%     At{i} = A.time;
%     Ad{i} = A.depth;
%     ASee{i} = A.See;
% end


%addpath('C:\Users\nsc4750\OneDrive - UNC-Wilmington\CMS Work\Fall (0825)\SmartMooring Data')
addpath('SmartMooring Data')
SMnames = {'roxsi_smartmooring_L2_E01sp_1851.mat','roxsi_smartmooring_L2_E02sp_1859.mat',...
    'roxsi_smartmooring_L2_E05sp_1853.mat','roxsi_smartmooring_L2_E07sp_1855.mat',...
    'roxsi_smartmooring_L2_E07sp_1857.mat','roxsi_smartmooring_L2_E08sp_1852.mat',...
    'roxsi_smartmooring_L2_E09sp_1850.mat','roxsi_smartmooring_L2_E09sp_1856.mat',...
    'roxsi_smartmooring_L2_E10sp_1848.mat','roxsi_smartmooring_L2_E11sp_1860.mat',...
    'roxsi_smartmooring_L2_E13sp_1849.mat'};
% for i = 1:length(SMnames)
%     load(SMnames{i})
%     SMlat{i} = spotsmartL2.latitude;
%     SMlon{i} = spotsmartL2.longitude;
%     SMt{i} = spotsmartL2.dtime;
%     SMd{i} = spotsmartL2.depth;
% 
% % .6059, .9632

%addpath('C:\Users\nsc4750\OneDrive - UNC-Wilmington\CMS Work\Summer (0516-0728)')
addpath('../Summer (0516-0728)/')
load('WBvariables.mat','Blat','Blon','Xlat','Xlon')


%% Plot of Asilomar Locations (only ADCP)

figure(1);clf;

    % Plot the Original Buoys (X01, X03, X04)
geoscatter(Xlat(1),Xlon(1),'r','o','LineWidth',5)
hold on
geoscatter(Xlat(2),Xlon(2),'ob','LineWidth',14)
geoscatter(Xlat(3),Xlon(3),'om','LineWidth',14)


    % Plot ADCPnames{3} (X05)
load(ADCPnames{3})
geoscatter(A.latitude,A.longitude,'xy','LineWidth',14)
title('Asilomar: Instrument Locations','fontsize',18)
    % Plot ADCPnames{4} (X11)
load(ADCPnames{4})
geoscatter(A.latitude,A.longitude,'xg','LineWidth',14)

geobasemap satellite
legend('X01 (h = 20.7m)','X03 (h = 15.8m)','X04 (h = 11.2m)',...
    'X05 (h = 9.3m)','X11 (h = 2.72m)','location','NorthEast','fontsize',12)


%% Plot of China Rock Locations (ADCP and smart spotters)

figure(2);clf;

    % Plot Original Buoys (B01, B03, B05)
geoscatter(Blat(1),Blon(1),'ro','LineWidth',10)
hold on
geoscatter(Blat(2),Blon(2),'bo','LineWidth',10)
geoscatter(Blat(3),Blon(3),'mo','LineWidth',10)

    % Plot ADCPnames{1} (B10)
load(ADCPnames{1})
geoscatter(A.latitude,A.longitude,'yx','LineWidth',10)
title('China Rock: Instrument Locations','fontsize',18)

    % Plot ADCPnames{2} (B13)
load(ADCPnames{2})
geoscatter(A.latitude,A.longitude,'gx','LineWidth',10)

geobasemap satellite

    % Plot SMnames{1} (E01)
% load(SMnames{1});
% geoscatter(spotsmartL2.latitude,spotsmartL2.longitude,'or','LineWidth',2)
%     % Plot SMnames{2} (E02)
% load(SMnames{2});
% geoscatter(spotsmartL2.latitude,spotsmartL2.longitude,'og','LineWidth',2)
%     % Plot SMnames{3} (E05)
% load(SMnames{3});
% geoscatter(spotsmartL2.latitude,spotsmartL2.longitude,'ob','LineWidth',2)
%     % Plot SMnames{4} and SMnames{5} (E07)
% load(SMnames{4});
% geoscatter(spotsmartL2.latitude,spotsmartL2.longitude,'oy','LineWidth',2)
%     % Plot SMnames{6} (E08)
% load(SMnames{6});
% geoscatter(spotsmartL2.latitude,spotsmartL2.longitude,'ok','LineWidth',2)
%     % Plot SMnames{7} and SMnames{8} (E09)
% load(SMnames{7});
% geoscatter(spotsmartL2.latitude,spotsmartL2.longitude,'om','LineWidth',2)
%     % Plot SMnames{9} (E10)
% load(SMnames{9});
% geoscatter(spotsmartL2.latitude,spotsmartL2.longitude,'oc','LineWidth',2)
%     % Plot SMnames{10} (E11)
% load(SMnames{10});
% geoscatter(spotsmartL2.latitude,spotsmartL2.longitude,'ow','LineWidth',2)
%     % Plot SMnames{11} (E13)
% load(SMnames{11});
% geoscatter(spotsmartL2.latitude,spotsmartL2.longitude,'or','LineWidth',2)


% legend('B10','B13','E01','E02','E05','E07','E08','E09','E10','E11',...
%     'E13','B01','B03','B05','location','northwest')

legend('B01 (h = 33.9m)','B03 (h = 21.0m)','B05 (h = 16.2m)',...
    'B10 (h = 9.5m)','B13 (h = 4.7m)','location','NorthEast','fontsize',12)


















