%% plot_TaylorDiags_All
%
% Noah Clark            7/25/24
%
%
% Purpose: Use the statistics comparing the 4 models (SWAN050, SWAN050sm,
%          COUP050, & SWAN075) to the observed measurements (calculated in
%          'calc_MvsO_Stats_All.m') to plot taylor diagrams. These taylor
%          diagrams represent the models over the majority of the frequency
%          spectrum as opposed to breaking them up into the 4 frequency
%          bands as done previously.
%           
%        - 6 diagrams total: Energy, Direction, Spread, TED, Difference
%          in Direction, & Difference in Spread


%% Preliminaries

clc;clear;

addpath('../taylor diagram/taylor diagram') % path for taylor diagram function
load('Stats_FS.mat')



%% Plot Formatting

% - STATS1: SWAN050
% - STATS2: SWAN050sm
% - STATS3: COUP050
% - STATS4: SWAN075
% - STATS: concatonated version of them all


% ----------------- FOR STATISTICS AT EACH BUOY -------------------------

    % Create the color list:
% - Color the markers based on which model
%           b = SWAN050
%           r = SWAN050sm
%           m = COUP050
%           g = SWAN075

colorList = {'b','b','b','b','b','b','b','b','b','r','r','r','r','r',...
    'r','r','r','r','m','m','m','m','m','m','m','m','m','g','g','g',...
    'g','g','g','g','g','g'};
FaceColorList = {'b','b','b','b','b','none','b','b','b','r','r','r','r','r',...
    'none','r','r','r','m','m','m','m','m','none','m','m','m','g','g','g',...
    'g','g','none','g','g','g'};

    % Create the marker symbols:
%           o = B01       o = X01
%           s = B03       x = X03
%           d = B05       + = X04
%           v = B10       * = X05
%           < = B13

MarkerType = {'o','s','d','v','<','o','x','+','*','o','s','d','v','<',...
    'o','x','+','*','o','s','d','v','<','o','x','+','*','o','s','d',...
    'v','<','o','x','+','*'};


    % Create the name list to label the legend:
NameList = {'SWAN050-B01','SWAN050-B03','SWAN050-B05','SWAN050-B10',...
    'SWAN050-B13','SWAN050-X01','SWAN050-X03','SWAN050-X04','SWAN050-X05',...
    'SWAN050sm-B01','SWAN050sm-B03','SWAN050sm-B05','SWAN050sm-B10',...
    'SWAN050sm-B13','SWAN050sm-X01','SWAN050sm-X03','SWAN050sm-X04',...
    'SWAN050sm-X05','COUP050-B01','COUP050-B03','COUP050-B05',...
    'COUP050-B10','COUP050-B13','COUP050-X01','COUP050-X03','COUP050-X04',...
    'COUP050-X05','SWAN075-B01','SWAN075-B03','SWAN075-B05','SWAN075-B10',...
    'SWAN075-B13','SWAN075-X01','SWAN075-X03','SWAN075-X04','SWAN075-X05'};



% ----------------- FOR STATISTICS BETWEEN EACH BUOY ----------------------

    % Create the color list:
% - Color the markers based on which model
%           b = SWAN050
%           r = SWAN050sm
%           m = COUP050
%           g = SWAN075

colorList_BT = {'b','b','b','b','b','b','b','r','r','r','r','r','r','r',...
    'm','m','m','m','m','m','m','g','g','g','g','g','g','g'};
FaceColorList_BT = {'b','b','b','b','none','b','b','r','r','r','r','none','r','r',...
    'm','m','m','m','none','m','m','g','g','g','g','none','g','g'};


    % Create the marker symbols:
%           o = B01-B03       o = X01-X03
%           s = B03-B05       x = X03-X04
%           d = B05-B10       + = X04-X05
%           v = B10-B13     

MarkerType_BT = {'o','s','d','v','o','x','+','o','s','d','v','o','x','+',...
    'o','s','d','v','o','x','+','o','s','d','v','o','x','+'};


    % Create the name list to label the legend:
NameList_BT = {'SWAN050:B01-B03','SWAN050:B03-B05','SWAN050:B05-B10',...
    'SWAN050:B10-B13','SWAN050:X01-X03','SWAN050:X03-X04','SWAN050:X04-X05',...
    'SWAN050sm:B01-B03','SWAN050sm:B03-B05','SWAN050sm:B05-B10',...
    'SWAN050sm:B10-B13','SWAN050sm:X01-X03','SWAN050sm:X03-X04','SWAN050sm:X04-X05',...
    'COUP050:B01-B03','COUP050:B03-B05','COUP050:B05-B10',...
    'COUP050:B10-B13','COUP050:X01-X03','COUP050:X03-X04','COUP050:X04-X05',...
    'SWAN075:B01-B03','SWAN075:B03-B05','SWAN075:B05-B10',...
    'SWAN075:B10-B13','SWAN075:X01-X03','SWAN075:X03-X04','SWAN075:X04-X05'};



%% Create Taylor Diagrams for the Statistics at Each Buoy
% 1x3 subplots (E, Dir, Spread)

    % Create Figure:
Taylor_FS = figure(1); clf;
set(gcf,'position',[150,443,1275,478])


% --------------------------- ENERGY --------------------------------------

STATS1 = [Stats_FS.EnergyStats.STD(1,:)./Stats_FS.iObsSTD.E_STD(1,:);...
    Stats_FS.EnergyStats.RMSE(1,:)./Stats_FS.iObsSTD.E_STD(1,:);...
    Stats_FS.EnergyStats.MaxCor(1,:)];
STATS2 = [Stats_FS.EnergyStats.STD(2,:)./Stats_FS.iObsSTD.E_STD(2,:);...
    Stats_FS.EnergyStats.RMSE(2,:)./Stats_FS.iObsSTD.E_STD(2,:);...
    Stats_FS.EnergyStats.MaxCor(2,:)];
STATS3 = [Stats_FS.EnergyStats.STD(3,:)./Stats_FS.iObsSTD.E_STD(3,:);...
    Stats_FS.EnergyStats.RMSE(3,:)./Stats_FS.iObsSTD.E_STD(3,:);...
    Stats_FS.EnergyStats.MaxCor(3,:)];
STATS4 = [Stats_FS.EnergyStats.STD(4,:)./Stats_FS.iObsSTD.E_STD(4,:);...
    Stats_FS.EnergyStats.RMSE(4,:)./Stats_FS.iObsSTD.E_STD(4,:);...
    Stats_FS.EnergyStats.MaxCor(4,:)];

STATS = cat(2,STATS1,STATS2,STATS3,STATS4);


subplot(1,3,1)
Tay1 = STaylorDiag(STATS);
for i = 1:length(STATS)
    Tay1.SPlot(STATS(1,i),STATS(2,i),STATS(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Energy (f: 0.04 - 0.23 Hz)',''})

clear STATS1 STATS2 STATS3 STATS4 STATS


% -------------------------- DIRECTION ------------------------------------

STATS1 = [Stats_FS.DirectionStats.STD(1,:)./Stats_FS.iObsSTD.Dir_STD(1,:);...
    Stats_FS.DirectionStats.RMSE(1,:)./Stats_FS.iObsSTD.Dir_STD(1,:);...
    Stats_FS.DirectionStats.MaxCor(1,:)];
STATS2 = [Stats_FS.DirectionStats.STD(2,:)./Stats_FS.iObsSTD.Dir_STD(2,:);...
    Stats_FS.DirectionStats.RMSE(2,:)./Stats_FS.iObsSTD.Dir_STD(2,:);...
    Stats_FS.DirectionStats.MaxCor(2,:)];
STATS3 = [Stats_FS.DirectionStats.STD(3,:)./Stats_FS.iObsSTD.Dir_STD(3,:);...
    Stats_FS.DirectionStats.RMSE(3,:)./Stats_FS.iObsSTD.Dir_STD(3,:);...
    Stats_FS.DirectionStats.MaxCor(3,:)];
STATS4 = [Stats_FS.DirectionStats.STD(4,:)./Stats_FS.iObsSTD.Dir_STD(4,:);...
    Stats_FS.DirectionStats.RMSE(4,:)./Stats_FS.iObsSTD.Dir_STD(4,:);...
    Stats_FS.DirectionStats.MaxCor(4,:)];

STATS = cat(2,STATS1,STATS2,STATS3,STATS4);


subplot(1,3,2)
Tay2 = STaylorDiag(STATS);
for i = 1:length(STATS)
    Tay2.SPlot(STATS(1,i),STATS(2,i),STATS(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Direction (f: 0.04 - 0.23 Hz)',''})

clear STATS1 STATS2 STATS3 STATS4 STATS

% --------------------------- SPREAD --------------------------------------

STATS1 = [Stats_FS.SpreadStats.STD(1,:)./Stats_FS.iObsSTD.Spread_STD(1,:);...
    Stats_FS.SpreadStats.RMSE(1,:)./Stats_FS.iObsSTD.Spread_STD(1,:);...
    Stats_FS.SpreadStats.MaxCor(1,:)];
STATS2 = [Stats_FS.SpreadStats.STD(2,:)./Stats_FS.iObsSTD.Spread_STD(2,:);...
    Stats_FS.SpreadStats.RMSE(2,:)./Stats_FS.iObsSTD.Spread_STD(2,:);...
    Stats_FS.SpreadStats.MaxCor(2,:)];
STATS3 = [Stats_FS.SpreadStats.STD(3,:)./Stats_FS.iObsSTD.Spread_STD(3,:);...
    Stats_FS.SpreadStats.RMSE(3,:)./Stats_FS.iObsSTD.Spread_STD(3,:);...
    Stats_FS.SpreadStats.MaxCor(3,:)];
STATS4 = [Stats_FS.SpreadStats.STD(4,:)./Stats_FS.iObsSTD.Spread_STD(4,:);...
    Stats_FS.SpreadStats.RMSE(4,:)./Stats_FS.iObsSTD.Spread_STD(4,:);...
    Stats_FS.SpreadStats.MaxCor(4,:)];

STATS = cat(2,STATS1,STATS2,STATS3,STATS4);


subplot(1,3,3)
Tay3 = STaylorDiag(STATS);
for i = 1:length(STATS)
    Tay3.SPlot(STATS(1,i),STATS(2,i),STATS(3,i),'Marker',MarkerType{i},...
        'MarkerSize',6,'Color',colorList{i},'MarkerFaceColor',FaceColorList{i});
end
title({'Spread (f: 0.04 - 0.23 Hz)',''})
legend(NameList,'FontSize',7,'NumColumns',9,'Location',...
    [0.091241840859396,0.031729428172943,0.83921567199277,0.120292883787195])

clear STATS1 STATS2 STATS3 STATS4 STATS


%% Create Taylor Diagrams for the Statistics BETWEEN Each Buoy Pair
% 1x3 subplots (TED, diffDir, diffSpread)

    % Create Figure:
Taylor_FS_BT = figure(2); clf;
set(gcf,'position',[300,443,1275,478])


% ------------------------------ TED --------------------------------------

STATS1 = [Stats_FS.TEDStats.STD(1,:)./Stats_FS.iObsSTD.TED_STD(1,:);...
    Stats_FS.TEDStats.RMSE(1,:)./Stats_FS.iObsSTD.TED_STD(1,:);...
    Stats_FS.TEDStats.MaxCor(1,:)];
STATS2 = [Stats_FS.TEDStats.STD(2,:)./Stats_FS.iObsSTD.TED_STD(2,:);...
    Stats_FS.TEDStats.RMSE(2,:)./Stats_FS.iObsSTD.TED_STD(2,:);...
    Stats_FS.TEDStats.MaxCor(2,:)];
STATS3 = [Stats_FS.TEDStats.STD(3,:)./Stats_FS.iObsSTD.TED_STD(3,:);...
    Stats_FS.TEDStats.RMSE(3,:)./Stats_FS.iObsSTD.TED_STD(3,:);...
    Stats_FS.TEDStats.MaxCor(3,:)];
STATS4 = [Stats_FS.TEDStats.STD(4,:)./Stats_FS.iObsSTD.TED_STD(4,:);...
    Stats_FS.TEDStats.RMSE(4,:)./Stats_FS.iObsSTD.TED_STD(4,:);...
    Stats_FS.TEDStats.MaxCor(4,:)];

STATS = cat(2,STATS1,STATS2,STATS3,STATS4);

subplot(1,3,1)
Tay1 = STaylorDiag(STATS);
for i = 1:length(STATS)
    Tay1.SPlot(STATS(1,i),STATS(2,i),STATS(3,i),'Marker',MarkerType_BT{i},...
        'MarkerSize',6,'Color',colorList_BT{i},'MarkerFaceColor',FaceColorList_BT{i});
end
title({'TED (f: 0.04 - 0.23 Hz)',''})

clear STATS1 STATS2 STATS3 STATS4 STATS

% -----------------------DIFFERENCE IN DIRECTION --------------------------

STATS1 = [Stats_FS.diffDirStats.STD(1,:)./Stats_FS.iObsSTD.diffDir_STD(1,:);...
    Stats_FS.diffDirStats.RMSE(1,:)./Stats_FS.iObsSTD.diffDir_STD(1,:);...
    Stats_FS.diffDirStats.MaxCor(1,:)];
STATS2 = [Stats_FS.diffDirStats.STD(2,:)./Stats_FS.iObsSTD.diffDir_STD(2,:);...
    Stats_FS.diffDirStats.RMSE(2,:)./Stats_FS.iObsSTD.diffDir_STD(2,:);...
    Stats_FS.diffDirStats.MaxCor(2,:)];
STATS3 = [Stats_FS.diffDirStats.STD(3,:)./Stats_FS.iObsSTD.diffDir_STD(3,:);...
    Stats_FS.diffDirStats.RMSE(3,:)./Stats_FS.iObsSTD.diffDir_STD(3,:);...
    Stats_FS.diffDirStats.MaxCor(3,:)];
STATS4 = [Stats_FS.diffDirStats.STD(4,:)./Stats_FS.iObsSTD.diffDir_STD(4,:);...
    Stats_FS.diffDirStats.RMSE(4,:)./Stats_FS.iObsSTD.diffDir_STD(4,:);...
    Stats_FS.diffDirStats.MaxCor(4,:)];

STATS = cat(2,STATS1,STATS2,STATS3,STATS4);

subplot(1,3,2)
Tay2 = STaylorDiag(STATS);
for i = 1:length(STATS)
    Tay2.SPlot(STATS(1,i),STATS(2,i),STATS(3,i),'Marker',MarkerType_BT{i},...
        'MarkerSize',6,'Color',colorList_BT{i},'MarkerFaceColor',FaceColorList_BT{i});
end
title({'Diff Direction (f: 0.04 - 0.23 Hz)',''})

clear STATS1 STATS2 STATS3 STATS4 STATS


% -----------------------DIFFERENCE IN SPREAD -----------------------------

STATS1 = [Stats_FS.diffSpreadStats.STD(1,:)./Stats_FS.iObsSTD.diffSpread_STD(1,:);...
    Stats_FS.diffSpreadStats.RMSE(1,:)./Stats_FS.iObsSTD.diffSpread_STD(1,:);...
    Stats_FS.diffSpreadStats.MaxCor(1,:)];
STATS2 = [Stats_FS.diffSpreadStats.STD(2,:)./Stats_FS.iObsSTD.diffSpread_STD(2,:);...
    Stats_FS.diffSpreadStats.RMSE(2,:)./Stats_FS.iObsSTD.diffSpread_STD(2,:);...
    Stats_FS.diffSpreadStats.MaxCor(2,:)];
STATS3 = [Stats_FS.diffSpreadStats.STD(3,:)./Stats_FS.iObsSTD.diffSpread_STD(3,:);...
    Stats_FS.diffSpreadStats.RMSE(3,:)./Stats_FS.iObsSTD.diffSpread_STD(3,:);...
    Stats_FS.diffSpreadStats.MaxCor(3,:)];
STATS4 = [Stats_FS.diffSpreadStats.STD(4,:)./Stats_FS.iObsSTD.diffSpread_STD(4,:);...
    Stats_FS.diffSpreadStats.RMSE(4,:)./Stats_FS.iObsSTD.diffSpread_STD(4,:);...
    Stats_FS.diffSpreadStats.MaxCor(4,:)];

STATS = cat(2,STATS1,STATS2,STATS3,STATS4);

subplot(1,3,3)
Tay3 = STaylorDiag(STATS);
for i = 1:length(STATS)
    Tay3.SPlot(STATS(1,i),STATS(2,i),STATS(3,i),'Marker',MarkerType_BT{i},...
        'MarkerSize',6,'Color',colorList_BT{i},'MarkerFaceColor',FaceColorList_BT{i});
end
title({'Diff Spread (f: 0.04 - 0.23 Hz)',''})

legend(NameList_BT,'FontSize',7,'NumColumns',9,'Location',...
    [0.091241840859396,0.031729428172943,0.83921567199277,0.120292883787195])

clear STATS1 STATS2 STATS3 STATS4 STATS



%% Save Figures & Clear Variables


saveas(Taylor_FS,'Figures/Taylor_FS.jpeg')
saveas(Taylor_FS_BT,'Figures/Taylor_FS_BT.jpeg')


clear Tay1 Tay2 Tay3 i















