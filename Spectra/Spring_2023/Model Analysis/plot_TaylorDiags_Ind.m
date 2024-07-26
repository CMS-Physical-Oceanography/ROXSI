%% plot_TaylorDiags_Ind.m
%
% Noah Clark            5/29/2024
%
%
% Purpose: Use the statistics comparing the model data to the observed
%          data, which was calculated in plot_MvsO.m and plot_MvsO_BT.m, 
%          to plot taylor diagrams of the statistics of each individual
%          buoy for each model.
%

% TODO:
%       - How do the ticks and labels look?



%% Preliminaries

clc;clear;

addpath('taylor diagram/taylor diagram')

load('obsVSmod_STATS.mat')
load('ModelResults.mat')
clear E_OSTD D_OSTD DS_OSTD


%% Plot Formatting

    % Create the color list:
% - Color the markers based on which model
%           b = SWAN050
%           r = SWAN050sm
%           m = COUP050
%           g = SWAN075

colorList = {'b','r','m','g','b','r','m','g','b','r','m','g','b','r',...
    'm','g','b','r','m','g','b','r','m','g','b','r','m','g','b','r',...
    'm','g','b','r','m','g'}; 
FaceColorList = {'b','r','m','g','b','r','m','g','b','r','m','g','b',...
    'r','m','g','b','r','m','g','none','none','none','none','b','r',...
    'm','g','b','r','m','g','b','r','m','g'};

    % Create the marker symbols:
%           o = B01       o = X01
%           s = B03       x = X03
%           d = B05       + = X04
%           v = B10       * = X05
%           < = B13

MarkerType = {'o','o','o','o','s','s','s','s','d','d','d','d','v','v','v','v',...
    '<','<','<','<','o','o','o','o','x','x','x','x','+','+','+','+','*','*','*','*'};

    % Create the name list to label the legend:
NameList = {'SWAN050-B01','SWAN050sm-B01','COUP050-B01','SWAN075-B01',...
    'SWAN050-B03','SWAN050sm-B03','COUP050-B03','SWAN075-B03',...
    'SWAN050-B05','SWAN050sm-B05','COUP050-B05','SWAN075-B05',...
    'SWAN050-B10','SWAN050sm-B10','COUP050-B10','SWAN075-B10',...
    'SWAN050-B13','SWAN050sm-B13','COUP050-B13','SWAN075-B13',...
    'SWAN050-X01','SWAN050sm-X01','COUP050-X01','SWAN075-X01',...
    'SWAN050-X03','SWAN050sm-X03','COUP050-X03','SWAN075-X03',...
    'SWAN050-X04','SWAN050sm-X04','COUP050-X04','SWAN075-X04',...
    'SWAN050-X05','SWAN050sm-X05','COUP050-X05','SWAN075-X05'};


%% Plot Energy Taylor Diagram

        % CREATE THE 'STATS' ARRAYS TO ENTER INTO TAYLOR DIAGRAM FUNCTION:
        
    % 1st Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 1;
        STATS1(:,i) = [EnergyStats.STD(M,a+1)/iObsSTD.EnergySTD(M,a), EnergyStats.RMSE(M,a+1)/iObsSTD.EnergySTD(M,a), EnergyStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 5;
        STATS1(:,i) = [EnergyStats.STD(M,b+1)/iObsSTD.EnergySTD(M,b), EnergyStats.RMSE(M,b+1)/iObsSTD.EnergySTD(M,b), EnergyStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 9;
        STATS1(:,i) = [EnergyStats.STD(M,c+1)/iObsSTD.EnergySTD(M,c), EnergyStats.RMSE(M,c+1)/iObsSTD.EnergySTD(M,c), EnergyStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 13;
        STATS1(:,i) = [EnergyStats.STD(M,d+1)/iObsSTD.EnergySTD(M,d), EnergyStats.RMSE(M,d+1)/iObsSTD.EnergySTD(M,d), EnergyStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 2nd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 2;
        STATS2(:,i) = [EnergyStats.STD(M,a+1)/iObsSTD.EnergySTD(M,a), EnergyStats.RMSE(M,a+1)/iObsSTD.EnergySTD(M,a), EnergyStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 6;
        STATS2(:,i) = [EnergyStats.STD(M,b+1)/iObsSTD.EnergySTD(M,b), EnergyStats.RMSE(M,b+1)/iObsSTD.EnergySTD(M,b), EnergyStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 10;
        STATS2(:,i) = [EnergyStats.STD(M,c+1)/iObsSTD.EnergySTD(M,c), EnergyStats.RMSE(M,c+1)/iObsSTD.EnergySTD(M,c), EnergyStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 14;
        STATS2(:,i) = [EnergyStats.STD(M,d+1)/iObsSTD.EnergySTD(M,d), EnergyStats.RMSE(M,d+1)/iObsSTD.EnergySTD(M,d), EnergyStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 3rd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 3;
        STATS3(:,i) = [EnergyStats.STD(M,a+1)/iObsSTD.EnergySTD(M,a), EnergyStats.RMSE(M,a+1)/iObsSTD.EnergySTD(M,a), EnergyStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 7;
        STATS3(:,i) = [EnergyStats.STD(M,b+1)/iObsSTD.EnergySTD(M,b), EnergyStats.RMSE(M,b+1)/iObsSTD.EnergySTD(M,b), EnergyStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 11;
        STATS3(:,i) = [EnergyStats.STD(M,c+1)/iObsSTD.EnergySTD(M,c), EnergyStats.RMSE(M,c+1)/iObsSTD.EnergySTD(M,c), EnergyStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 15;
        STATS3(:,i) = [EnergyStats.STD(M,d+1)/iObsSTD.EnergySTD(M,d), EnergyStats.RMSE(M,d+1)/iObsSTD.EnergySTD(M,d), EnergyStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 4th Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 4;
        STATS4(:,i) = [EnergyStats.STD(M,a+1)/iObsSTD.EnergySTD(M,a), EnergyStats.RMSE(M,a+1)/iObsSTD.EnergySTD(M,a), EnergyStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 8;
        STATS4(:,i) = [EnergyStats.STD(M,b+1)/iObsSTD.EnergySTD(M,b), EnergyStats.RMSE(M,b+1)/iObsSTD.EnergySTD(M,b), EnergyStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 12;
        STATS4(:,i) = [EnergyStats.STD(M,c+1)/iObsSTD.EnergySTD(M,c), EnergyStats.RMSE(M,c+1)/iObsSTD.EnergySTD(M,c), EnergyStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 16;
        STATS4(:,i) = [EnergyStats.STD(M,d+1)/iObsSTD.EnergySTD(M,d), EnergyStats.RMSE(M,d+1)/iObsSTD.EnergySTD(M,d), EnergyStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

                                % PLOT:
Taylor_i = figure(1);
clf;
set(gcf,'position',[144,48,1184,946])

subplot(3,4,1) % 1st frequency band
Tay1 = STaylorDiag(STATS1);
for i = 1:length(STATS1)
    Tay1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Energy: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,2)
Tay2 = STaylorDiag(STATS2);
for i = 1:length(STATS2)
    Tay2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Energy: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,3)
Tay3 = STaylorDiag(STATS3);
for i = 1:length(STATS3)
    Tay3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Energy: (0.8>f/f_p\geq1.2)',''})

subplot(3,4,4)
Tay4 = STaylorDiag(STATS4);
for i = 1:length(STATS4)
    Tay4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Energy: (1.2>f/f_p\geq2.0)',''})
legend(NameList,'FontSize',7,'location','north','NumColumns',9,...
    'location',[0.019646199700776,0.020483643891815,0.96537160440474,0.060782239376616])


iEnergy_STATS{1} = STATS1;
iEnergy_STATS{2} = STATS2;
iEnergy_STATS{3} = STATS3;
iEnergy_STATS{4} = STATS4;

clear STATS1 STATS2 STATS3 STATS4



%% Plot Direction Taylor Diagram

        % CREATE THE 'STATS' ARRAYS TO ENTER INTO TAYLOR DIAGRAM FUNCTION:

    % 1st Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 1;
        STATS1(:,i) = [DirectionStats.STD(M,a+1)/iObsSTD.DirSTD(M,a), DirectionStats.RMSE(M,a+1)/iObsSTD.DirSTD(M,a), DirectionStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 5;
        STATS1(:,i) = [DirectionStats.STD(M,b+1)/iObsSTD.DirSTD(M,b), DirectionStats.RMSE(M,b+1)/iObsSTD.DirSTD(M,b), DirectionStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 9;
        STATS1(:,i) = [DirectionStats.STD(M,c+1)/iObsSTD.DirSTD(M,c), DirectionStats.RMSE(M,c+1)/iObsSTD.DirSTD(M,c), DirectionStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 13;
        STATS1(:,i) = [DirectionStats.STD(M,d+1)/iObsSTD.DirSTD(M,d), DirectionStats.RMSE(M,d+1)/iObsSTD.DirSTD(M,d), DirectionStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 2nd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 2;
        STATS2(:,i) = [DirectionStats.STD(M,a+1)/iObsSTD.DirSTD(M,a), DirectionStats.RMSE(M,a+1)/iObsSTD.DirSTD(M,a), DirectionStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 6;
        STATS2(:,i) = [DirectionStats.STD(M,b+1)/iObsSTD.DirSTD(M,b), DirectionStats.RMSE(M,b+1)/iObsSTD.DirSTD(M,b), DirectionStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 10;
        STATS2(:,i) = [DirectionStats.STD(M,c+1)/iObsSTD.DirSTD(M,c), DirectionStats.RMSE(M,c+1)/iObsSTD.DirSTD(M,c), DirectionStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 14;
        STATS2(:,i) = [DirectionStats.STD(M,d+1)/iObsSTD.DirSTD(M,d), DirectionStats.RMSE(M,d+1)/iObsSTD.DirSTD(M,d), DirectionStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 3rd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 3;
        STATS3(:,i) = [DirectionStats.STD(M,a+1)/iObsSTD.DirSTD(M,a), DirectionStats.RMSE(M,a+1)/iObsSTD.DirSTD(M,a), DirectionStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 7;
        STATS3(:,i) = [DirectionStats.STD(M,b+1)/iObsSTD.DirSTD(M,b), DirectionStats.RMSE(M,b+1)/iObsSTD.DirSTD(M,b), DirectionStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 11;
        STATS3(:,i) = [DirectionStats.STD(M,c+1)/iObsSTD.DirSTD(M,c), DirectionStats.RMSE(M,c+1)/iObsSTD.DirSTD(M,c), DirectionStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 15;
        STATS3(:,i) = [DirectionStats.STD(M,d+1)/iObsSTD.DirSTD(M,d), DirectionStats.RMSE(M,d+1)/iObsSTD.DirSTD(M,d), DirectionStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 4th Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 4;
        STATS4(:,i) = [DirectionStats.STD(M,a+1)/iObsSTD.DirSTD(M,a), DirectionStats.RMSE(M,a+1)/iObsSTD.DirSTD(M,a), DirectionStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 8;
        STATS4(:,i) = [DirectionStats.STD(M,b+1)/iObsSTD.DirSTD(M,b), DirectionStats.RMSE(M,b+1)/iObsSTD.DirSTD(M,b), DirectionStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 12;
        STATS4(:,i) = [DirectionStats.STD(M,c+1)/iObsSTD.DirSTD(M,c), DirectionStats.RMSE(M,c+1)/iObsSTD.DirSTD(M,c), DirectionStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 16;
        STATS4(:,i) = [DirectionStats.STD(M,d+1)/iObsSTD.DirSTD(M,d), DirectionStats.RMSE(M,d+1)/iObsSTD.DirSTD(M,d), DirectionStats.MaxCor(M,d+1)];
        d = d+1;
    end
end


                                % PLOT:
%Taylor_iDir = figure(2);
%clf;
%set(gcf,'position',[532,252,1121,503])

subplot(3,4,5) % 1st frequency band
Tay1 = STaylorDiag(STATS1);
for i = 1:length(STATS1)
    Tay1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Direction: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,6)
Tay2 = STaylorDiag(STATS2);
for i = 1:length(STATS2)
    Tay2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Direction: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,7)
Tay3 = STaylorDiag(STATS3);
for i = 1:length(STATS3)
    Tay3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Direction: (0.8>f/f_p\geq1.2)',''})

subplot(3,4,8)
Tay4 = STaylorDiag(STATS4);
for i = 1:length(STATS4)
    Tay4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Direction: (1.2>f/f_p\geq2.0)',''})
%legend(NameList,'FontSize',7,'location','north','NumColumns',9,...
%    'location',[0.106733192084153,0.245592813823111,0.811467889908257,0.058766859344892])


iDirection_STATS{1} = STATS1;
iDirection_STATS{2} = STATS2;
iDirection_STATS{3} = STATS3;
iDirection_STATS{4} = STATS4;

clear STATS1 STATS2 STATS3 STATS4


%% Plot Directional Spread Taylor Diagram

        % CREATE THE 'STATS' ARRAYS TO ENTER INTO TAYLOR DIAGRAM FUNCTION:

    % 1st Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 1;
        STATS1(:,i) = [DSpreadStats.STD(M,a+1)/iObsSTD.DSpreadSTD(M,a), DSpreadStats.RMSE(M,a+1)/iObsSTD.DSpreadSTD(M,a), DSpreadStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 5;
        STATS1(:,i) = [DSpreadStats.STD(M,b+1)/iObsSTD.DSpreadSTD(M,b), DSpreadStats.RMSE(M,b+1)/iObsSTD.DSpreadSTD(M,b), DSpreadStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 9;
        STATS1(:,i) = [DSpreadStats.STD(M,c+1)/iObsSTD.DSpreadSTD(M,c), DSpreadStats.RMSE(M,c+1)/iObsSTD.DSpreadSTD(M,c), DSpreadStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 13;
        STATS1(:,i) = [DirectionStats.STD(M,d+1)/iObsSTD.DSpreadSTD(M,d), DirectionStats.RMSE(M,d+1)/iObsSTD.DSpreadSTD(M,d), DirectionStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 2nd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 2;
        STATS2(:,i) = [DSpreadStats.STD(M,a+1)/iObsSTD.DSpreadSTD(M,a), DSpreadStats.RMSE(M,a+1)/iObsSTD.DSpreadSTD(M,a), DSpreadStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 6;
        STATS2(:,i) = [DSpreadStats.STD(M,b+1)/iObsSTD.DSpreadSTD(M,b), DSpreadStats.RMSE(M,b+1)/iObsSTD.DSpreadSTD(M,b), DSpreadStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 10;
        STATS2(:,i) = [DSpreadStats.STD(M,c+1)/iObsSTD.DSpreadSTD(M,c), DSpreadStats.RMSE(M,c+1)/iObsSTD.DSpreadSTD(M,c), DSpreadStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 14;
        STATS2(:,i) = [DirectionStats.STD(M,d+1)/iObsSTD.DSpreadSTD(M,d), DirectionStats.RMSE(M,d+1)/iObsSTD.DSpreadSTD(M,d), DirectionStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 3rd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 3;
        STATS3(:,i) = [DSpreadStats.STD(M,a+1)/iObsSTD.DSpreadSTD(M,a), DSpreadStats.RMSE(M,a+1)/iObsSTD.DSpreadSTD(M,a), DSpreadStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 7;
        STATS3(:,i) = [DSpreadStats.STD(M,b+1)/iObsSTD.DSpreadSTD(M,b), DSpreadStats.RMSE(M,b+1)/iObsSTD.DSpreadSTD(M,b), DSpreadStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 11;
        STATS3(:,i) = [DSpreadStats.STD(M,c+1)/iObsSTD.DSpreadSTD(M,c), DSpreadStats.RMSE(M,c+1)/iObsSTD.DSpreadSTD(M,c), DSpreadStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 15;
        STATS3(:,i) = [DirectionStats.STD(M,d+1)/iObsSTD.DSpreadSTD(M,d), DirectionStats.RMSE(M,d+1)/iObsSTD.DSpreadSTD(M,d), DirectionStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 4th Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:36 %(4models*9buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 || i==29 || i==33
        M = 4;
        STATS4(:,i) = [DSpreadStats.STD(M,a+1)/iObsSTD.DSpreadSTD(M,a), DSpreadStats.RMSE(M,a+1)/iObsSTD.DSpreadSTD(M,a), DSpreadStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 || i==30 || i==34
        M = 8;
        STATS4(:,i) = [DSpreadStats.STD(M,b+1)/iObsSTD.DSpreadSTD(M,b), DSpreadStats.RMSE(M,b+1)/iObsSTD.DSpreadSTD(M,b), DSpreadStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 || i==31 || i==35
        M = 12;
        STATS4(:,i) = [DSpreadStats.STD(M,c+1)/iObsSTD.DSpreadSTD(M,c), DSpreadStats.RMSE(M,c+1)/iObsSTD.DSpreadSTD(M,c), DSpreadStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 16;
        STATS4(:,i) = [DirectionStats.STD(M,d+1)/iObsSTD.DSpreadSTD(M,d), DirectionStats.RMSE(M,d+1)/iObsSTD.DSpreadSTD(M,d), DirectionStats.MaxCor(M,d+1)];
        d = d+1;
    end
end


                                % PLOT:
% Taylor_iDSpread = figure(3);
% clf;
% set(gcf,'position',[532,32,1121,503])

subplot(3,4,9) % 1st frequency band
Tay1 = STaylorDiag(STATS1);
for i = 1:length(STATS1)
    Tay1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'D Spread: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,10)
Tay2 = STaylorDiag(STATS2);
for i = 1:length(STATS2)
    Tay2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'D Spread: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,11)
Tay3 = STaylorDiag(STATS3);
for i = 1:length(STATS3)
    Tay3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'D Spread: (0.8>f/f_p\geq1.2)',''})

subplot(3,4,12)
Tay4 = STaylorDiag(STATS4);
for i = 1:length(STATS4)
    Tay4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'D Spread: (1.2>f/f_p\geq2.0)',''})
%legend(NameList,'FontSize',7,'location','north','NumColumns',9,...
%    'location',[0.106733192084153,0.245592813823111,0.811467889908257,0.058766859344892])


iDSpread_STATS{1} = STATS1;
iDSpread_STATS{2} = STATS2;
iDSpread_STATS{3} = STATS3;
iDSpread_STATS{4} = STATS4;

clear STATS1 STATS2 STATS3 STATS4




%% ----------------------------------------------------------------
%%              BETWEEN BUOY STATISTICS
%% ----------------------------------------------------------------



%% Plot Formatting (for b/t buoy stats)

    % Create the color list:
% - Color the markers based on which model
%           b = SWAN050
%           r = SWAN050sm
%           m = COUP050
%           g = SWAN075

colorList = {'b','r','m','g','b','r','m','g','b','r','m','g','b','r',...
    'm','g','b','r','m','g','b','r','m','g','b','r','m','g'}; 
FaceColorList = {'b','r','m','g','b','r','m','g','b','r','m','g','b','r',...
    'm','g','none','none','none','none','b','r','m','g','b','r','m','g'};

    % Create the marker symbols:
%           o = B01-B03       o = X01-X03
%           s = B03-B05       x = X03-X04
%           d = B05-B10       + = X04-X05
%           v = B10-B13       

MarkerType = {'o','o','o','o','s','s','s','s','d','d','d','d','v','v','v','v',...
    'o','o','o','o','x','x','x','x','+','+','+','+'};

    % Create the name list to label the legend:
NameList = {'SWAN050:B01-B03','SWAN050sm:B01-B03','COUP050:B01-B03','SWAN075:B01-B03',...
    'SWAN050:B03-B04','SWAN050sm:B03-B04','COUP050:B03-B04','SWAN075:B03-B04',...
    'SWAN050:B05-B10','SWAN050sm:B05-B10','COUP050:B05-B10','SWAN075:B05-B10',...
    'SWAN050:B10-B13','SWAN050sm:B10-B13','COUP050:B10-B13','SWAN075:B10-B13',...
    'SWAN050:X01-X03','SWAN050sm:X01-X03','COUP050:X01-X03','SWAN075:X01-X03',...
    'SWAN050:X03-X04','SWAN050sm:X03-X04','COUP050:X03-X04','SWAN075:X03-X04',...
    'SWAN050:X04-X05','SWAN050sm:X04-X05','COUP050:X04-X05','SWAN075:X04-X05'};




%% Plot TED Taylor Diagram

        % CREATE THE 'STATS' ARRAYS TO ENTER INTO TAYLOR DIAGRAM FUNCTION:
        
    % 1st Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25
        M = 1;
        STATS1(:,i) = [TEDStats.STD(M,a+1)/iObsSTD_BT.TEDi_STD(M,a), TEDStats.RMSE(M,a+1)/iObsSTD_BT.TEDi_STD(M,a), TEDStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26
        M = 5;
        STATS1(:,i) = [TEDStats.STD(M,b+1)/iObsSTD_BT.TEDi_STD(M,b), TEDStats.RMSE(M,b+1)/iObsSTD_BT.TEDi_STD(M,b), TEDStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27
        M = 9;
        STATS1(:,i) = [TEDStats.STD(M,c+1)/iObsSTD_BT.TEDi_STD(M,c), TEDStats.RMSE(M,c+1)/iObsSTD_BT.TEDi_STD(M,c), TEDStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 13;
        STATS1(:,i) = [TEDStats.STD(M,d+1)/iObsSTD_BT.TEDi_STD(M,d), TEDStats.RMSE(M,d+1)/iObsSTD_BT.TEDi_STD(M,d), TEDStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 2nd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25
        M = 2;
        STATS2(:,i) = [TEDStats.STD(M,a+1)/iObsSTD_BT.TEDi_STD(M,a), TEDStats.RMSE(M,a+1)/iObsSTD_BT.TEDi_STD(M,a), TEDStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26
        M = 6;
        STATS2(:,i) = [TEDStats.STD(M,b+1)/iObsSTD_BT.TEDi_STD(M,b), TEDStats.RMSE(M,b+1)/iObsSTD_BT.TEDi_STD(M,b), TEDStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 
        M = 10;
        STATS2(:,i) = [TEDStats.STD(M,c+1)/iObsSTD_BT.TEDi_STD(M,c), TEDStats.RMSE(M,c+1)/iObsSTD_BT.TEDi_STD(M,c), TEDStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 14;
        STATS2(:,i) = [TEDStats.STD(M,d+1)/iObsSTD_BT.TEDi_STD(M,d), TEDStats.RMSE(M,d+1)/iObsSTD_BT.TEDi_STD(M,d), TEDStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 3rd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 
        M = 3;
        STATS3(:,i) = [TEDStats.STD(M,a+1)/iObsSTD_BT.TEDi_STD(M,a), TEDStats.RMSE(M,a+1)/iObsSTD_BT.TEDi_STD(M,a), TEDStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26
        M = 7;
        STATS3(:,i) = [TEDStats.STD(M,b+1)/iObsSTD_BT.TEDi_STD(M,b), TEDStats.RMSE(M,b+1)/iObsSTD_BT.TEDi_STD(M,b), TEDStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27
        M = 11;
        STATS3(:,i) = [TEDStats.STD(M,c+1)/iObsSTD_BT.TEDi_STD(M,c), TEDStats.RMSE(M,c+1)/iObsSTD_BT.TEDi_STD(M,c), TEDStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 15;
        STATS3(:,i) = [TEDStats.STD(M,d+1)/iObsSTD_BT.TEDi_STD(M,d), TEDStats.RMSE(M,d+1)/iObsSTD_BT.TEDi_STD(M,d), TEDStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 4th Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 
        M = 4;
        STATS4(:,i) = [TEDStats.STD(M,a+1)/iObsSTD_BT.TEDi_STD(M,a), TEDStats.RMSE(M,a+1)/iObsSTD_BT.TEDi_STD(M,a), TEDStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 
        M = 8;
        STATS4(:,i) = [TEDStats.STD(M,b+1)/iObsSTD_BT.TEDi_STD(M,b), TEDStats.RMSE(M,b+1)/iObsSTD_BT.TEDi_STD(M,b), TEDStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27
        M = 12;
        STATS4(:,i) = [TEDStats.STD(M,c+1)/iObsSTD_BT.TEDi_STD(M,c), TEDStats.RMSE(M,c+1)/iObsSTD_BT.TEDi_STD(M,c), TEDStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 16;
        STATS4(:,i) = [TEDStats.STD(M,d+1)/iObsSTD_BT.TEDi_STD(M,d), TEDStats.RMSE(M,d+1)/iObsSTD_BT.TEDi_STD(M,d), TEDStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

                                % PLOT:
Taylor_i_BT = figure(2);
clf;
set(gcf,'position',[350,48,1184,946])

subplot(3,4,1) % 1st frequency band
Tay1 = STaylorDiag(STATS1);
for i = 1:length(STATS1)
    Tay1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'TED: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,2)
Tay2 = STaylorDiag(STATS2);
for i = 1:length(STATS2)
    Tay2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'TED: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,3)
Tay3 = STaylorDiag(STATS3);
for i = 1:length(STATS3)
    Tay3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'TED: (0.8>f/f_p\geq1.2)',''})

subplot(3,4,4)
Tay4 = STaylorDiag(STATS4);
for i = 1:length(STATS4)
    Tay4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'TED: (1.2>f/f_p\geq2.0)',''})
legend(NameList,'FontSize',7,'location','north','NumColumns',9,...
    'location',[0.019646199700776,0.020483643891815,0.96537160440474,0.060782239376616])


iTED_STATS{1} = STATS1;
iTED_STATS{2} = STATS2;
iTED_STATS{3} = STATS3;
iTED_STATS{4} = STATS4;

clear STATS1 STATS2 STATS3 STATS4



%% Plot Difference in Direction B/T Buoys Taylor Diags


        % CREATE THE 'STATS' ARRAYS TO ENTER INTO TAYLOR DIAGRAM FUNCTION:
        
    % 1st Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25
        M = 1;
        STATS1(:,i) = [diffDirStats.STD(M,a+1)/iObsSTD_BT.dDi_STD(M,a), diffDirStats.RMSE(M,a+1)/iObsSTD_BT.dDi_STD(M,a), diffDirStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26
        M = 5;
        STATS1(:,i) = [diffDirStats.STD(M,b+1)/iObsSTD_BT.dDi_STD(M,b), diffDirStats.RMSE(M,b+1)/iObsSTD_BT.dDi_STD(M,b), diffDirStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27
        M = 9;
        STATS1(:,i) = [diffDirStats.STD(M,c+1)/iObsSTD_BT.dDi_STD(M,c), diffDirStats.RMSE(M,c+1)/iObsSTD_BT.dDi_STD(M,c), diffDirStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 13;
        STATS1(:,i) = [diffDirStats.STD(M,d+1)/iObsSTD_BT.dDi_STD(M,d), diffDirStats.RMSE(M,d+1)/iObsSTD_BT.dDi_STD(M,d), diffDirStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 2nd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25
        M = 2;
        STATS2(:,i) = [diffDirStats.STD(M,a+1)/iObsSTD_BT.dDi_STD(M,a), diffDirStats.RMSE(M,a+1)/iObsSTD_BT.dDi_STD(M,a), diffDirStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26
        M = 6;
        STATS2(:,i) = [diffDirStats.STD(M,b+1)/iObsSTD_BT.dDi_STD(M,b), diffDirStats.RMSE(M,b+1)/iObsSTD_BT.dDi_STD(M,b), diffDirStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 
        M = 10;
        STATS2(:,i) = [diffDirStats.STD(M,c+1)/iObsSTD_BT.dDi_STD(M,c), diffDirStats.RMSE(M,c+1)/iObsSTD_BT.dDi_STD(M,c), diffDirStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 14;
        STATS2(:,i) = [diffDirStats.STD(M,d+1)/iObsSTD_BT.dDi_STD(M,d), diffDirStats.RMSE(M,d+1)/iObsSTD_BT.dDi_STD(M,d), diffDirStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 3rd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 
        M = 3;
        STATS3(:,i) = [diffDirStats.STD(M,a+1)/iObsSTD_BT.dDi_STD(M,a), diffDirStats.RMSE(M,a+1)/iObsSTD_BT.dDi_STD(M,a), diffDirStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26
        M = 7;
        STATS3(:,i) = [diffDirStats.STD(M,b+1)/iObsSTD_BT.dDi_STD(M,b), diffDirStats.RMSE(M,b+1)/iObsSTD_BT.dDi_STD(M,b), diffDirStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27
        M = 11;
        STATS3(:,i) = [diffDirStats.STD(M,c+1)/iObsSTD_BT.dDi_STD(M,c), diffDirStats.RMSE(M,c+1)/iObsSTD_BT.dDi_STD(M,c), diffDirStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 15;
        STATS3(:,i) = [diffDirStats.STD(M,d+1)/iObsSTD_BT.dDi_STD(M,d), diffDirStats.RMSE(M,d+1)/iObsSTD_BT.dDi_STD(M,d), diffDirStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 4th Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 
        M = 4;
        STATS4(:,i) = [diffDirStats.STD(M,a+1)/iObsSTD_BT.dDi_STD(M,a), diffDirStats.RMSE(M,a+1)/iObsSTD_BT.dDi_STD(M,a), diffDirStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 
        M = 8;
        STATS4(:,i) = [diffDirStats.STD(M,b+1)/iObsSTD_BT.dDi_STD(M,b), diffDirStats.RMSE(M,b+1)/iObsSTD_BT.dDi_STD(M,b), diffDirStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27
        M = 12;
        STATS4(:,i) = [diffDirStats.STD(M,c+1)/iObsSTD_BT.dDi_STD(M,c), diffDirStats.RMSE(M,c+1)/iObsSTD_BT.dDi_STD(M,c), diffDirStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 16;
        STATS4(:,i) = [diffDirStats.STD(M,d+1)/iObsSTD_BT.dDi_STD(M,d), diffDirStats.RMSE(M,d+1)/iObsSTD_BT.dDi_STD(M,d), diffDirStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

                                % PLOT:

subplot(3,4,5) % 1st frequency band
Tay1 = STaylorDiag(STATS1);
for i = 1:length(STATS1)
    Tay1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Diff Dir: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,6)
Tay2 = STaylorDiag(STATS2);
for i = 1:length(STATS2)
    Tay2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Diff Dir: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,7)
Tay3 = STaylorDiag(STATS3);
for i = 1:length(STATS3)
    Tay3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Diff Dir: (0.8>f/f_p\geq1.2)',''})

subplot(3,4,8)
Tay4 = STaylorDiag(STATS4);
for i = 1:length(STATS4)
    Tay4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Diff Dir: (1.2>f/f_p\geq2.0)',''})


idiffD_STATS{1} = STATS1;
idiffD_STATS{2} = STATS2;
idiffD_STATS{3} = STATS3;
idiffD_STATS{4} = STATS4;

clear STATS1 STATS2 STATS3 STATS4



%% Plot Difference in Spread B/T Buoys Taylor Diags


        % CREATE THE 'STATS' ARRAYS TO ENTER INTO TAYLOR DIAGRAM FUNCTION:
        
    % 1st Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25
        M = 1;
        STATS1(:,i) = [diffDspreadStats.STD(M,a+1)/iObsSTD_BT.dDSi_STD(M,a), diffDspreadStats.RMSE(M,a+1)/iObsSTD_BT.dDSi_STD(M,a), diffDspreadStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26
        M = 5;
        STATS1(:,i) = [diffDspreadStats.STD(M,b+1)/iObsSTD_BT.dDSi_STD(M,b), diffDspreadStats.RMSE(M,b+1)/iObsSTD_BT.dDSi_STD(M,b), diffDspreadStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27
        M = 9;
        STATS1(:,i) = [diffDspreadStats.STD(M,c+1)/iObsSTD_BT.dDSi_STD(M,c), diffDspreadStats.RMSE(M,c+1)/iObsSTD_BT.dDSi_STD(M,c), diffDspreadStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 13;
        STATS1(:,i) = [diffDspreadStats.STD(M,d+1)/iObsSTD_BT.dDSi_STD(M,d), diffDspreadStats.RMSE(M,d+1)/iObsSTD_BT.dDSi_STD(M,d), diffDspreadStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 2nd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25
        M = 2;
        STATS2(:,i) = [diffDspreadStats.STD(M,a+1)/iObsSTD_BT.dDSi_STD(M,a), diffDspreadStats.RMSE(M,a+1)/iObsSTD_BT.dDSi_STD(M,a), diffDspreadStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26
        M = 6;
        STATS2(:,i) = [diffDspreadStats.STD(M,b+1)/iObsSTD_BT.dDSi_STD(M,b), diffDspreadStats.RMSE(M,b+1)/iObsSTD_BT.dDSi_STD(M,b), diffDspreadStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27 
        M = 10;
        STATS2(:,i) = [diffDspreadStats.STD(M,c+1)/iObsSTD_BT.dDSi_STD(M,c), diffDspreadStats.RMSE(M,c+1)/iObsSTD_BT.dDSi_STD(M,c), diffDspreadStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 14;
        STATS2(:,i) = [diffDspreadStats.STD(M,d+1)/iObsSTD_BT.dDSi_STD(M,d), diffDspreadStats.RMSE(M,d+1)/iObsSTD_BT.dDSi_STD(M,d), diffDspreadStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 3rd Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 
        M = 3;
        STATS3(:,i) = [diffDspreadStats.STD(M,a+1)/iObsSTD_BT.dDSi_STD(M,a), diffDspreadStats.RMSE(M,a+1)/iObsSTD_BT.dDSi_STD(M,a), diffDspreadStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26
        M = 7;
        STATS3(:,i) = [diffDspreadStats.STD(M,b+1)/iObsSTD_BT.dDSi_STD(M,b), diffDspreadStats.RMSE(M,b+1)/iObsSTD_BT.dDSi_STD(M,b), diffDspreadStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27
        M = 11;
        STATS3(:,i) = [diffDspreadStats.STD(M,c+1)/iObsSTD_BT.dDSi_STD(M,c), diffDspreadStats.RMSE(M,c+1)/iObsSTD_BT.dDSi_STD(M,c), diffDspreadStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 15;
        STATS3(:,i) = [diffDspreadStats.STD(M,d+1)/iObsSTD_BT.dDSi_STD(M,d), diffDspreadStats.RMSE(M,d+1)/iObsSTD_BT.dDSi_STD(M,d), diffDspreadStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

    % 4th Frequency Band
a = 1; b = 1; c = 1; d = 1;
for i = 1:28 %(4models*7 b/t buoys)
    if i==1 || i==5 || i==9 || i==13|| i==17 || i==21 || i==25 
        M = 4;
        STATS4(:,i) = [diffDspreadStats.STD(M,a+1)/iObsSTD_BT.dDSi_STD(M,a), diffDspreadStats.RMSE(M,a+1)/iObsSTD_BT.dDSi_STD(M,a), diffDspreadStats.MaxCor(M,a+1)];
        a = a+1;
    elseif i==2 || i==6 || i==10 || i==14 || i==18 || i==22 || i==26 
        M = 8;
        STATS4(:,i) = [diffDspreadStats.STD(M,b+1)/iObsSTD_BT.dDSi_STD(M,b), diffDspreadStats.RMSE(M,b+1)/iObsSTD_BT.dDSi_STD(M,b), diffDspreadStats.MaxCor(M,b+1)];
        b = b+1;
    elseif i == 3 || i==7 || i==11 || i==15 || i==19 || i==23 || i==27
        M = 12;
        STATS4(:,i) = [diffDspreadStats.STD(M,c+1)/iObsSTD_BT.dDSi_STD(M,c), diffDspreadStats.RMSE(M,c+1)/iObsSTD_BT.dDSi_STD(M,c), diffDspreadStats.MaxCor(M,c+1)];
        c = c+1;
    else
        M = 16;
        STATS4(:,i) = [diffDspreadStats.STD(M,d+1)/iObsSTD_BT.dDSi_STD(M,d), diffDspreadStats.RMSE(M,d+1)/iObsSTD_BT.dDSi_STD(M,d), diffDspreadStats.MaxCor(M,d+1)];
        d = d+1;
    end
end

                                % PLOT:

subplot(3,4,9) % 1st frequency band
Tay1 = STaylorDiag(STATS1);
for i = 1:length(STATS1)
    Tay1.SPlot(STATS1(1,i),STATS1(2,i),STATS1(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Diff Spread: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,10)
Tay2 = STaylorDiag(STATS2);
for i = 1:length(STATS2)
    Tay2.SPlot(STATS2(1,i),STATS2(2,i),STATS2(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Diff Spread: (0.5>f/f_p\geq0.8)',''})

subplot(3,4,11)
Tay3 = STaylorDiag(STATS3);
for i = 1:length(STATS3)
    Tay3.SPlot(STATS3(1,i),STATS3(2,i),STATS3(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Diff Spread: (0.8>f/f_p\geq1.2)',''})

subplot(3,4,12)
Tay4 = STaylorDiag(STATS4);
for i = 1:length(STATS4)
    Tay4.SPlot(STATS4(1,i),STATS4(2,i),STATS4(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Diff Spread: (1.2>f/f_p\geq2.0)',''})


idiffDS_STATS{1} = STATS1;
idiffDS_STATS{2} = STATS2;
idiffDS_STATS{3} = STATS3;
idiffDS_STATS{4} = STATS4;

clear STATS1 STATS2 STATS3 STATS4



%% ----------------------------------------------------------------
%%              Save Figures & Clear Variables
%% ----------------------------------------------------------------


saveas(figure(1),'Figures/TaylorDiags/Taylor_i.jpeg')
saveas(figure(2),'Figures/TaylorDiags/Taylor_i_BT.jpeg')



clear a ans b c d i M Tay1 Tay2 Tay3 Tay4















