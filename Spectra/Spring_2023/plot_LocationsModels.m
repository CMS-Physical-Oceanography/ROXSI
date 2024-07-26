%% plot_LocationsModels.m
%
% Noah Clark        4/25/2024
%
% Purpose: To plot the locations of the spotter buoys and ADCPs at Asilomar
%          and China Rock. Also, to plot the locations of the model's 
%          observation points.


clc;clear;

%% Load in Data

    % Load ADCP Data
addpath('../Fall (0825)/ADCP Data')
ADCPnames = {'roxsi_signature_L2_B10_See.mat',...
    'roxsi_signature_L2_B13_See.mat','roxsi_signature_L2_X05_See.mat',...
    'roxsi_signature_L2_X11_See.mat'};

    % Load Spotter Results
addpath('../')
load('WBvariables.mat','Blat','Blon','Xlat','Xlon')


%% Load Model Results, Determine Time-Averages , and Re-save
    % Load Model Results
addpath('Model Results')
load('SWAN2022_Rogkn050_sub.mat')

B01.AVGspec = mean(B01.WS,2);
B03.AVGspec = mean(B03.WS,2);
B05.AVGspec = mean(B05.WS,2);
B10.AVGspec = mean(B10.WS,2);
B13.AVGspec = mean(B13.WS,2);
X01.AVGspec = mean(X01.WS,2);
X03.AVGspec = mean(X03.WS,2);
X04.AVGspec = mean(X04.WS,2);
X05.AVGspec = mean(X05.WS,2);

% Add Results into structure and resave a new mat file
Model.B01 = B01;
Model.B03 = B03;
Model.B05 = B05;
Model.B07 = B07;
Model.B10 = B10;
Model.B13 = B13;
Model.X01 = X01;
Model.X03 = X03;
Model.X04 = X04;
Model.X05 = X05;
clear B01 B03 B05 B07 B10 B13 X01 X03 X04 X05 X08
% save('ModelResults.mat','Model')


%% Plot of Asilomar Locations (only ADCP)

figure(1);clf;

    % Plot the Original Buoys (X01, X03, X04)
geoscatter(Xlat(1),Xlon(1),'r','o','LineWidth',10)
hold on
geoscatter(Xlat(2),Xlon(2),'ob','LineWidth',10)
geoscatter(Xlat(3),Xlon(3),'om','LineWidth',10)

    % Plot X05 (ADCPnames{3}) 
load(ADCPnames{3})
geoscatter(A.latitude,A.longitude,'xy','LineWidth',14)
title('Asilomar: Instrument Locations','fontsize',18)
    % Plot X11 (ADCPnames{4})
load(ADCPnames{4})
geoscatter(A.latitude,A.longitude,'xg','LineWidth',14)

    % Plot Model Locations
geoscatter(Model.X01.Lat,Model.X01.Lon,'^r','LineWidth',3)
geoscatter(Model.X03.Lat,Model.X03.Lon,'^b','LineWidth',3)
geoscatter(Model.X04.Lat,Model.X04.Lon,'^m','LineWidth',3)
geoscatter(Model.X05.Lat,Model.X05.Lon,'^y','LineWidth',3)
geoscatter(Model.X08.Lat,Model.X08.Lon,'^g','LineWidth',3)

geobasemap satellite

legend('X01 (h = 20.7m)','X03 (h = 15.8m)','X04 (h = 11.2m)',...
    'X05 (h = 9.3m)','X11 (h = 2.72m)','Model X01','Model X03',...
    'Model X04','Model X05','Model X08','location','NorthEast',...
    'NumColumns',2,'fontsize',12)

ax = gca;
ax.FontSize = 12;



%% Plot of China Rock Locations (ADCP and smart spotters)

figure(2);clf;

    % Plot Original Buoys (B01, B03, B05)
geoscatter(Blat(1),Blon(1),'ro','LineWidth',10)
hold on
geoscatter(Blat(2),Blon(2),'bo','LineWidth',10)
geoscatter(Blat(3),Blon(3),'mo','LineWidth',10)

    % Plot B10 (ADCPnames{1})
load(ADCPnames{1})
geoscatter(A.latitude,A.longitude,'yx','LineWidth',14)
title('China Rock: Instrument Locations','fontsize',20)

    % Plot B13 (ADCPnames{2})
load(ADCPnames{2})
geoscatter(A.latitude,A.longitude,'gx','LineWidth',14)

    % Plot Model Locations
geoscatter(Model.B01.Lat,Model.B01.Lon,'^r','LineWidth',3)
geoscatter(Model.B03.Lat,Model.B03.Lon,'^b','LineWidth',3)
geoscatter(Model.B05.Lat,Model.B05.Lon,'^m','LineWidth',3)
geoscatter(Model.B10.Lat,Model.B10.Lon,'^y','LineWidth',3)
geoscatter(Model.B13.Lat,Model.B13.Lon,'^g','LineWidth',3)


geobasemap satellite

legend('B01 (h = 33.9m)','B03 (h = 21.0m)','B05 (h = 16.2m)',...
    'B10 (h = 9.5m)','B13 (h = 4.7m)','Model B01','Model B03',...
    'Model B05','Model B10','Model B13','location','NorthEast',...
    'NumColumns',2,'fontsize',12)

ax = gca;
ax.FontSize = 12;



% THIS DOES NOTHING (WHY?)
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [7 7]);














