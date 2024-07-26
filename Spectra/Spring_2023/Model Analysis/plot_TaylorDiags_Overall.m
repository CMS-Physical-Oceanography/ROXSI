%% plot_TaylorDiags_Overall.m
%
% Noah Clark            5/29/2024
%
% Purpose: Use the statistics comparing the model data to the observed
%          data, which was calculated in plot_MvsO.m, to plot taylor
%          diagrams of the statistis of all buoys combined.
%
% Script Walkthrough: 
%   #1) Create the statistics arrays containing the energy statistics.
%       Input the statistic arrays into the STaylorDiag function to create
%       Taylor diagrams for the 4 different frequency ranges.
%   #2) Exact same as #1 but with direction.
%   #3) Exact same as #1 but with directional spread.
%   #4) Exact same as #1 but with TED.
%   #5) Exact same as #1 but with difference in direction between buoys.
%   #6) Exact same as #1 but with difference in directional spread between buoys.
%
% NOTE: In the future we plan to combine all of the subplots for each
%       figure into one plot. We may even try to combine all of the
%       plots in each figure into one master Taylor diagram.


%% Preliminaries

clc;clear;

addpath('taylor diagram/taylor diagram')

load('obsVSmod_STATS.mat')
load('ModelResults.mat')

%% Note
% To normalize by the std and put the observed points on top of each other

% STATS2(:,1) = [E_OSTD{2}/E_OSTD{2} 0 1]; % Observed E stats for 2nd frequency band
% STATS2(:,2) = [EnergyStats.STD(2)/E_OSTD{2}, EnergyStats.RMSE(2)/E_OSTD{2}, EnergyStats.MaxCor(2)]; % SWAN E Stats for 2nd freq band
% STATS2(:,3) = [EnergyStats.STD(6)/E_OSTD{2}, EnergyStats.RMSE(6)/E_OSTD{2}, EnergyStats.MaxCor(6)]; % SWANsm E Stats for 2nd freq band
% STATS2(:,4) = [EnergyStats.STD(10)/E_OSTD{2}, EnergyStats.RMSE(10)/E_OSTD{2}, EnergyStats.MaxCor(10)]; % COUP E Stats for 2nd freq band

%% Plot Formatting

colorList = [0.3569    0.0784    0.0784
    0.6784    0.4471    0.1725
    0.1020    0.3882    0.5176
    0.9725    0.4196    0.4392];
MarkerType = {'o','diamond','pentagram','^','v'};

%% #1) Energy Statistic Figure

    % Make the arrays of statistics:
% order to enter: [std rmse CC]
STATS1(:,1) = [E_OSTD{1} 0 1]; % Observed E stats for 1st frequency band
STATS1(:,2) = [EnergyStats.STD(1,1), EnergyStats.RMSE(1,1), EnergyStats.MaxCor(1,1)]; % SWAN E Stats for 1st freq band
STATS1(:,3) = [EnergyStats.STD(5,1), EnergyStats.RMSE(5,1), EnergyStats.MaxCor(5,1)]; % SWANsm E Stats for 1st freq band
STATS1(:,4) = [EnergyStats.STD(9,1), EnergyStats.RMSE(9,1), EnergyStats.MaxCor(9,1)]; % COUP E Stats for 1st freq band

STATS2(:,1) = [E_OSTD{2} 0 1]; % Observed E stats for 2nd frequency band
STATS2(:,2) = [EnergyStats.STD(2,1), EnergyStats.RMSE(2,1), EnergyStats.MaxCor(2,1)]; % SWAN E Stats for 2nd freq band
STATS2(:,3) = [EnergyStats.STD(6,1), EnergyStats.RMSE(6,1), EnergyStats.MaxCor(6,1)]; % SWANsm E Stats for 2nd freq band
STATS2(:,4) = [EnergyStats.STD(10,1), EnergyStats.RMSE(10,1), EnergyStats.MaxCor(10,1)]; % COUP E Stats for 2nd freq band

STATS3(:,1) = [E_OSTD{3} 0 1]; % Observed E stats for 3rd frequency band
STATS3(:,2) = [EnergyStats.STD(3,1), EnergyStats.RMSE(3,1), EnergyStats.MaxCor(3,1)]; % SWAN E Stats for 3rd freq band
STATS3(:,3) = [EnergyStats.STD(7,1), EnergyStats.RMSE(7,1), EnergyStats.MaxCor(7,1)]; % SWANsm E Stats for 3rd freq band
STATS3(:,4) = [EnergyStats.STD(11,1), EnergyStats.RMSE(11,1), EnergyStats.MaxCor(11,1)]; % COUP E Stats for 3rd freq band

STATS4(:,1) = [E_OSTD{4} 0 1]; % Observed E stats for 4th frequency band
STATS4(:,2) = [EnergyStats.STD(4,1), EnergyStats.RMSE(4,1), EnergyStats.MaxCor(4,1)]; % SWAN E Stats for 4th freq band
STATS4(:,3) = [EnergyStats.STD(8,1), EnergyStats.RMSE(8,1), EnergyStats.MaxCor(8,1)]; % SWANsm E Stats for 4th freq band
STATS4(:,4) = [EnergyStats.STD(12,1), EnergyStats.RMSE(12,1), EnergyStats.MaxCor(12,1)]; % COUP E Stats for 4th freq band

Energy_STATS{1} = STATS1;
Energy_STATS{2} = STATS2;
Energy_STATS{3} = STATS3;
Energy_STATS{4} = STATS4;

    % Plot:
Taylor_Energy = figure('Units','normalized','Position',[0.0292 0.3056 0.7802 0.5435]);
clf;

subplot(1,4,1) % 1st frequency band
TD1 = STaylorDiag(STATS1);
for i = 1:size(STATS1,2)
    TD1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Energy: (0.5>f/f_p\geq0.8)',''})

subplot(1,4,2) % 2nd frequency band
TD2 = STaylorDiag(STATS2);
for i = 1:size(STATS2,2)
    TD2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Energy: (0.8>f/f_p\geq1.2)',''})

subplot(1,4,3) % 3rd frequency band
TD3 = STaylorDiag(STATS3);
for i = 1:size(STATS3,2)
    TD3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Energy: (1.2>f/f_p\geq2.0)',''})

subplot(1,4,4) % 4th frequency band
TD4 = STaylorDiag(STATS4);
for i = 1:size(STATS4,2)
    TD4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Energy: (2.0>f/f_p\geq3.0)',''})


clear STATS1 STATS2 STATS3 STATS4


%% #2) Direction Statistics Figure

    % Make the arrays of statistics:
% order to enter: [std rmse CC]
STATS1(:,1) = [D_OSTD{1} 0 1]; % Observed Dir stats for 1st frequency band
STATS1(:,2) = [DirectionStats.STD(1,1), DirectionStats.RMSE(1,1), DirectionStats.MaxCor(1,1)]; % SWAN Dir Stats for 1st freq band
STATS1(:,3) = [DirectionStats.STD(5,1), DirectionStats.RMSE(5,1), DirectionStats.MaxCor(5,1)]; % SWANsm Dir Stats for 1st freq band
STATS1(:,4) = [DirectionStats.STD(9,1), DirectionStats.RMSE(9,1), DirectionStats.MaxCor(9,1)]; % COUP Dir Stats for 1st freq band

STATS2(:,1) = [D_OSTD{2} 0 1]; % Observed Dir stats for 2nd frequency band
STATS2(:,2) = [DirectionStats.STD(2,1), DirectionStats.RMSE(2,1), DirectionStats.MaxCor(2,1)]; % SWAN Dir Stats for 2nd freq band
STATS2(:,3) = [DirectionStats.STD(6,1), DirectionStats.RMSE(6,1), DirectionStats.MaxCor(6,1)]; % SWANsm Dir Stats for 2nd freq band
STATS2(:,4) = [DirectionStats.STD(10,1), DirectionStats.RMSE(10,1), DirectionStats.MaxCor(10,1)]; % COUP Dir Stats for 2nd freq band

STATS3(:,1) = [D_OSTD{3} 0 1]; % Observed Dir stats for 3rd frequency band
STATS3(:,2) = [DirectionStats.STD(3,1), DirectionStats.RMSE(3,1), DirectionStats.MaxCor(3,1)]; % SWAN Dir Stats for 3rd freq band
STATS3(:,3) = [DirectionStats.STD(7,1), DirectionStats.RMSE(7,1), DirectionStats.MaxCor(7,1)]; % SWANsm Dir Stats for 3rd freq band
STATS3(:,4) = [DirectionStats.STD(11,1), DirectionStats.RMSE(11,1), DirectionStats.MaxCor(11,1)]; % COUP Dir Stats for 3rd freq band

STATS4(:,1) = [D_OSTD{4} 0 1]; % Observed Dir stats for 4th frequency band
STATS4(:,2) = [DirectionStats.STD(4,1), DirectionStats.RMSE(4,1), DirectionStats.MaxCor(4,1)]; % SWAN Dir Stats for 4th freq band
STATS4(:,3) = [DirectionStats.STD(8,1), DirectionStats.RMSE(8,1), DirectionStats.MaxCor(8,1)]; % SWANsm Dir Stats for 4th freq band
STATS4(:,4) = [DirectionStats.STD(12,1), DirectionStats.RMSE(12,1), DirectionStats.MaxCor(12,1)]; % COUP Dir Stats for 4th freq band

Direction_STATS{1} = STATS1;
Direction_STATS{2} = STATS2;
Direction_STATS{3} = STATS3;
Direction_STATS{4} = STATS4;


    % Plot:
Taylor_Dir = figure('Units','normalized','Position',[0.0492 0.3056 0.7802 0.5435]);
clf;

subplot(1,4,1) % 1st frequency band
TD1 = STaylorDiag(STATS1);
for i = 1:size(STATS1,2)
    TD1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Direction: (0.5>f/f_p\geq0.8)',''})

subplot(1,4,2) % 2nd frequency band
TD2 = STaylorDiag(STATS2);
for i = 1:size(STATS2,2)
    TD2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Direction: (0.8>f/f_p\geq1.2)',''})

subplot(1,4,3) % 3rd frequency band
TD3 = STaylorDiag(STATS3);
for i = 1:size(STATS3,2)
    TD3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Direction: (1.2>f/f_p\geq2.0)',''})

subplot(1,4,4) % 4th frequency band
TD4 = STaylorDiag(STATS4);
for i = 1:size(STATS4,2)
    TD4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Direction: (2.0>f/f_p\geq3.0)',''})


clear STATS1 STATS2 STATS3 STATS4


%% #3) Directional Spread Statistics Figure

    % Make the arrays of statistics:
% order to enter: [std rmse CC]
STATS1(:,1) = [DS_OSTD{1} 0 1]; % Observed Spread stats for 1st frequency band
STATS1(:,2) = [DSpreadStats.STD(1,1), DSpreadStats.RMSE(1,1), DSpreadStats.MaxCor(1,1)]; % SWAN Spread Stats for 1st freq band
STATS1(:,3) = [DSpreadStats.STD(5,1), DSpreadStats.RMSE(5,1), DSpreadStats.MaxCor(5,1)]; % SWANsm Spread Stats for 1st freq band
STATS1(:,4) = [DSpreadStats.STD(9,1), DSpreadStats.RMSE(9,1), DSpreadStats.MaxCor(9,1)]; % COUP Spread Stats for 1st freq band

STATS2(:,1) = [DS_OSTD{2} 0 1]; % Observed Spread stats for 2nd frequency band
STATS2(:,2) = [DSpreadStats.STD(2,1), DSpreadStats.RMSE(2,1), DSpreadStats.MaxCor(2,1)]; % SWAN Spread Stats for 2nd freq band
STATS2(:,3) = [DSpreadStats.STD(6,1), DSpreadStats.RMSE(6,1), DSpreadStats.MaxCor(6,1)]; % SWANsm Spread Stats for 2nd freq band
STATS2(:,4) = [DSpreadStats.STD(10,1), DSpreadStats.RMSE(10,1), DSpreadStats.MaxCor(10,1)]; % COUP Spread Stats for 2nd freq band

STATS3(:,1) = [DS_OSTD{3} 0 1]; % Observed Spread stats for 3rd frequency band
STATS3(:,2) = [DSpreadStats.STD(3,1), DSpreadStats.RMSE(3,1), DSpreadStats.MaxCor(3,1)]; % SWAN Spread Stats for 3rd freq band
STATS3(:,3) = [DSpreadStats.STD(7,1), DSpreadStats.RMSE(7,1), DSpreadStats.MaxCor(7,1)]; % SWANsm Spread Stats for 3rd freq band
STATS3(:,4) = [DSpreadStats.STD(11,1), DSpreadStats.RMSE(11,1), DSpreadStats.MaxCor(11,1)]; % COUP Spread Stats for 3rd freq band

STATS4(:,1) = [DS_OSTD{4} 0 1]; % Observed Spread stats for 4th frequency band
STATS4(:,2) = [DSpreadStats.STD(4,1), DSpreadStats.RMSE(4,1), DSpreadStats.MaxCor(4,1)]; % SWAN Spread Stats for 4th freq band
STATS4(:,3) = [DSpreadStats.STD(8,1), DSpreadStats.RMSE(8,1), DSpreadStats.MaxCor(8,1)]; % SWANsm Spread Stats for 4th freq band
STATS4(:,4) = [DSpreadStats.STD(12,1), DSpreadStats.RMSE(12,1), DSpreadStats.MaxCor(12,1)]; % COUP Spread Stats for 4th freq band

DSpread_STATS{1} = STATS1;
DSpread_STATS{2} = STATS2;
DSpread_STATS{3} = STATS3;
DSpread_STATS{4} = STATS4;


    % Plot:
Taylor_Dspread = figure('Units','normalized','Position',[0.0692 0.3056 0.7802 0.5435]);
clf;

subplot(1,4,1) % 1st frequency band
TD1 = STaylorDiag(STATS1);
for i = 1:size(STATS1,2)
    TD1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Spread: (0.5>f/f_p\geq0.8)',''})

subplot(1,4,2) % 2nd frequency band
TD2 = STaylorDiag(STATS2);
for i = 1:size(STATS2,2)
    TD2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Spread: (0.8>f/f_p\geq1.2)',''})

subplot(1,4,3) % 3rd frequency band
TD3 = STaylorDiag(STATS3);
for i = 1:size(STATS3,2)
    TD3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Spread: (1.2>f/f_p\geq2.0)',''})

subplot(1,4,4) % 4th frequency band
TD4 = STaylorDiag(STATS4);
for i = 1:size(STATS4,2)
    TD4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'Spread: (2.0>f/f_p\geq3.0)',''})


clear STATS1 STATS2 STATS3 STATS4


%% #4) TED Statistics Figure

% NOTE: - SWAN050sm has negative correlation for 1st freq band so it's off of plot
%       - COUP050 has negative correlation for 3rd freq band so it's off of plot

    % Make the arrays of statistics:
% order to enter: [std rmse CC]
STATS1(:,1) = [TED_OSTD{1} 0 1]; % Observed TED stats for 1st frequency band
STATS1(:,2) = [TEDStats.STD(1), TEDStats.RMSE(1), TEDStats.MaxCor(1)]; % SWAN Dir Stats for 1st freq band
STATS1(:,3) = [TEDStats.STD(5), TEDStats.RMSE(5), TEDStats.MaxCor(5)]; % SWANsm Dir Stats for 1st freq band
STATS1(:,4) = [TEDStats.STD(9), TEDStats.RMSE(9), TEDStats.MaxCor(9)]; % COUP Dir Stats for 1st freq band

STATS2(:,1) = [TED_OSTD{2} 0 1]; % Observed TED stats for 2nd frequency band
STATS2(:,2) = [TEDStats.STD(2), TEDStats.RMSE(2), TEDStats.MaxCor(2)]; % SWAN Dir Stats for 2nd freq band
STATS2(:,3) = [TEDStats.STD(6), TEDStats.RMSE(6), TEDStats.MaxCor(6)]; % SWANsm Dir Stats for 2nd freq band
STATS2(:,4) = [TEDStats.STD(10), TEDStats.RMSE(10), TEDStats.MaxCor(10)]; % COUP Dir Stats for 2nd freq band

STATS3(:,1) = [TED_OSTD{3} 0 1]; % Observed TED stats for 3rd frequency band
STATS3(:,2) = [TEDStats.STD(3), TEDStats.RMSE(3), TEDStats.MaxCor(3)]; % SWAN Dir Stats for 3rd freq band
STATS3(:,3) = [TEDStats.STD(7), TEDStats.RMSE(7), TEDStats.MaxCor(7)]; % SWANsm Dir Stats for 3rd freq band
STATS3(:,4) = [TEDStats.STD(11), TEDStats.RMSE(11), TEDStats.MaxCor(11)]; % COUP Dir Stats for 3rd freq band

STATS4(:,1) = [TED_OSTD{4} 0 1]; % Observed TED stats for 4th frequency band
STATS4(:,2) = [TEDStats.STD(4), TEDStats.RMSE(4), TEDStats.MaxCor(4)]; % SWAN Dir Stats for 4th freq band
STATS4(:,3) = [TEDStats.STD(8), TEDStats.RMSE(8), TEDStats.MaxCor(8)]; % SWANsm Dir Stats for 4th freq band
STATS4(:,4) = [TEDStats.STD(12), TEDStats.RMSE(12), TEDStats.MaxCor(12)]; % COUP Dir Stats for 4th freq band

TED_STATS{1} = STATS1;
TED_STATS{2} = STATS2;
TED_STATS{3} = STATS3;
TED_STATS{4} = STATS4;


    % Plot:
Taylor_TED = figure('Units','normalized','Position',[0.0892 0.3056 0.7802 0.5435]);

subplot(1,4,1) % 1st frequency band
TD1 = STaylorDiag(STATS1);
for i = 1:size(STATS1,2)
    TD1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'TED: (0.5>f/f_p\geq0.8)',''})

subplot(1,4,2) % 2nd frequency band
TD2 = STaylorDiag(STATS2);
for i = 1:size(STATS2,2)
    TD2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'TED: (0.8>f/f_p\geq1.2)',''})

subplot(1,4,3) % 3rd frequency band
TD3 = STaylorDiag(STATS3);
for i = 1:size(STATS3,2)
    TD3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'TED: (1.2>f/f_p\geq2.0)',''})

subplot(1,4,4) % 4th frequency band
TD4 = STaylorDiag(STATS4);
for i = 1:size(STATS4,2)
    TD4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'TED: (2.0>f/f_p\geq3.0)',''})


clear STATS1 STATS2 STATS3 STATS4



%% #5) Difference in Direction Statistics Figure

    % Make the arrays of statistics:
% order to enter: [std rmse CC]
STATS1(:,1) = [diffDir_OSTD{1} 0 1]; % Observed diffDir stats for 1st frequency band
STATS1(:,2) = [diffDirStats.STD(1), diffDirStats.RMSE(1), diffDirStats.MaxCor(1)]; % SWAN diffDir Stats for 1st freq band
STATS1(:,3) = [diffDirStats.STD(5), diffDirStats.RMSE(5), diffDirStats.MaxCor(5)]; % SWANsm diffDir Stats for 1st freq band
STATS1(:,4) = [diffDirStats.STD(9), diffDirStats.RMSE(9), diffDirStats.MaxCor(9)]; % COUP diffDir Stats for 1st freq band

STATS2(:,1) = [diffDir_OSTD{2} 0 1]; % Observed diffDir stats for 2nd frequency band
STATS2(:,2) = [diffDirStats.STD(2), diffDirStats.RMSE(2), diffDirStats.MaxCor(2)]; % SWAN diffDir Stats for 2nd freq band
STATS2(:,3) = [diffDirStats.STD(6), diffDirStats.RMSE(6), diffDirStats.MaxCor(6)]; % SWANsm diffDir Stats for 2nd freq band
STATS2(:,4) = [diffDirStats.STD(10), diffDirStats.RMSE(10), diffDirStats.MaxCor(10)]; % COUP diffDir Stats for 2nd freq band

STATS3(:,1) = [diffDir_OSTD{3} 0 1]; % Observed diffDir stats for 3rd frequency band
STATS3(:,2) = [diffDirStats.STD(3), diffDirStats.RMSE(3), diffDirStats.MaxCor(3)]; % SWAN diffDir Stats for 3rd freq band
STATS3(:,3) = [diffDirStats.STD(7), diffDirStats.RMSE(7), diffDirStats.MaxCor(7)]; % SWANsm diffDir Stats for 3rd freq band
STATS3(:,4) = [diffDirStats.STD(11), diffDirStats.RMSE(11), diffDirStats.MaxCor(11)]; % COUP diffDir Stats for 3rd freq band

STATS4(:,1) = [diffDir_OSTD{4} 0 1]; % Observed diffDir stats for 4th frequency band
STATS4(:,2) = [diffDirStats.STD(4), diffDirStats.RMSE(4), diffDirStats.MaxCor(4)]; % SWAN diffDir Stats for 4th freq band
STATS4(:,3) = [diffDirStats.STD(8), diffDirStats.RMSE(8), diffDirStats.MaxCor(8)]; % SWANsm diffDir Stats for 4th freq band
STATS4(:,4) = [diffDirStats.STD(12), diffDirStats.RMSE(12), diffDirStats.MaxCor(12)]; % COUP diffDir Stats for 4th freq band

diffDir_STATS{1} = STATS1;
diffDir_STATS{2} = STATS2;
diffDir_STATS{3} = STATS3;
diffDir_STATS{4} = STATS4;


    % Plot:
Taylor_diffDir = figure('Units','normalized','Position',[0.1092 0.3056 0.7802 0.5435]);

subplot(1,4,1) % 1st frequency band
TD1 = STaylorDiag(STATS1);
for i = 1:size(STATS1,2)
    TD1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'diffDir: (0.5>f/f_p\geq0.8)',''})

subplot(1,4,2) % 2nd frequency band
TD2 = STaylorDiag(STATS2);
for i = 1:size(STATS2,2)
    TD2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'diffDir: (0.8>f/f_p\geq1.2)',''})

subplot(1,4,3) % 3rd frequency band
TD3 = STaylorDiag(STATS3);
for i = 1:size(STATS3,2)
    TD3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'diffDir: (1.2>f/f_p\geq2.0)',''})

subplot(1,4,4) % 4th frequency band
TD4 = STaylorDiag(STATS4);
for i = 1:size(STATS4,2)
    TD4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'diffDir: (2.0>f/f_p\geq3.0)',''})


clear STATS1 STATS2 STATS3 STATS4



%% #6) Difference in Spread Statistics Figure


    % Make the arrays of statistics:
% order to enter: [std rmse CC]
STATS1(:,1) = [diffDspread_OSTD{1} 0 1]; % Observed diffDspread stats for 1st frequency band
STATS1(:,2) = [diffDspreadStats.STD(1), diffDspreadStats.RMSE(1), diffDspreadStats.MaxCor(1)]; % SWAN diffDspread Stats for 1st freq band
STATS1(:,3) = [diffDspreadStats.STD(5), diffDspreadStats.RMSE(5), diffDspreadStats.MaxCor(5)]; % SWANsm diffDspread Stats for 1st freq band
STATS1(:,4) = [diffDspreadStats.STD(9), diffDspreadStats.RMSE(9), diffDspreadStats.MaxCor(9)]; % COUP diffDspread Stats for 1st freq band

STATS2(:,1) = [diffDspread_OSTD{2} 0 1]; % Observed diffDspread stats for 2nd frequency band
STATS2(:,2) = [diffDspreadStats.STD(2), diffDspreadStats.RMSE(2), diffDspreadStats.MaxCor(2)]; % SWAN diffDspread Stats for 2nd freq band
STATS2(:,3) = [diffDspreadStats.STD(6), diffDspreadStats.RMSE(6), diffDspreadStats.MaxCor(6)]; % SWANsm diffDspread Stats for 2nd freq band
STATS2(:,4) = [diffDspreadStats.STD(10), diffDspreadStats.RMSE(10), diffDspreadStats.MaxCor(10)]; % COUP diffDspread Stats for 2nd freq band

STATS3(:,1) = [diffDspread_OSTD{3} 0 1]; % Observed diffDspread stats for 3rd frequency band
STATS3(:,2) = [diffDspreadStats.STD(3), diffDspreadStats.RMSE(3), diffDspreadStats.MaxCor(3)]; % SWAN diffDspread Stats for 3rd freq band
STATS3(:,3) = [diffDspreadStats.STD(7), diffDspreadStats.RMSE(7), diffDspreadStats.MaxCor(7)]; % SWANsm diffDspread Stats for 3rd freq band
STATS3(:,4) = [diffDspreadStats.STD(11), diffDspreadStats.RMSE(11), diffDspreadStats.MaxCor(11)]; % COUP diffDspread Stats for 3rd freq band

STATS4(:,1) = [diffDspread_OSTD{4} 0 1]; % Observed diffDspread stats for 4th frequency band
STATS4(:,2) = [diffDspreadStats.STD(4), diffDspreadStats.RMSE(4), diffDspreadStats.MaxCor(4)]; % SWAN diffDspread Stats for 4th freq band
STATS4(:,3) = [diffDspreadStats.STD(8), diffDspreadStats.RMSE(8), diffDspreadStats.MaxCor(8)]; % SWANsm diffDspread Stats for 4th freq band
STATS4(:,4) = [diffDspreadStats.STD(12), diffDspreadStats.RMSE(12), diffDspreadStats.MaxCor(12)]; % COUP diffDspread Stats for 4th freq band

diffDspread_STATS{1} = STATS1;
diffDspread_STATS{2} = STATS2;
diffDspread_STATS{3} = STATS3;
diffDspread_STATS{4} = STATS4;


    % Plot:
Taylor_diffDspread = figure('Units','normalized','Position',[0.1292 0.3056 0.7802 0.5435]);

subplot(1,4,1) % 1st frequency band
TD1 = STaylorDiag(STATS1);
for i = 1:size(STATS1,2)
    TD1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'diffDspread: (0.5>f/f_p\geq0.8)',''})

subplot(1,4,2) % 2nd frequency band
TD2 = STaylorDiag(STATS2);
for i = 1:size(STATS2,2)
    TD2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'diffDspread: (0.8>f/f_p\geq1.2)',''})

subplot(1,4,3) % 3rd frequency band
TD3 = STaylorDiag(STATS3);
for i = 1:size(STATS3,2)
    TD3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'diffDspread: (1.2>f/f_p\geq2.0)',''})

subplot(1,4,4) % 4th frequency band
TD4 = STaylorDiag(STATS4);
for i = 1:size(STATS4,2)
    TD4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},'MarkerSize',15,...
        'Color',colorList(i,:),'MarkerFaceColor',colorList(i,:));  
end
NameList = {'Observed','SWAN050','SWAN050sm','COUP050'};
legend(NameList,'FontSize',7,'location','best')
title({'diffDspread: (2.0>f/f_p\geq3.0)',''})


clear STATS1 STATS2 STATS3 STATS4





%% Save Plots

% saveas(Taylor_Energy,'Figures/TaylorDiags/Taylor_Energy.jpeg')
% saveas(Taylor_Dir,'Figures/TaylorDiags/Taylor_Dir.jpeg')
% saveas(Taylor_Dspread,'Figures/TaylorDiags/Taylor_Dspread.jpeg')
% saveas(Taylor_TED,'Figures/TaylorDiags/Taylor_TED.jpeg')
% saveas(Taylor_diffDir,'Figures/TaylorDiags/Taylor_diffDir.jpeg')
% saveas(Taylor_diffDspread,'Figures/TaylorDiags/Taylor_diffDspread.jpeg')












