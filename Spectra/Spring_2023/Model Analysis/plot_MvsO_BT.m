%% plot_MvsO_BT.m
%
% Noah Clark        6/3/2024
%
% Purpose: To determine the energy weighted mean TED, difference in
%          direction between buoys, and difference in directional spread
%          between the buoys for the models and for the observed in each of
%          the 4 frequency bands. Then plot each of the modeled directions
%          and spreads against the observed.
%
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%
% TO DO :
%       - finalize plots
%       - calculate all statistics for each one of buoys and make a third
%         dimension of the tables to save the individual stats in
%       - add the new model
%       - the last 4 aren't and haven't been working for ALL


%% Preliminaries

clc;clear;

load('InterpModelObs.mat')

addpath('../Stat Functions') % the functions for the statistical analysis
% the function for TED analysis:
%addpath('/Users/noahclark/Library/CloudStorage/OneDrive-UNC-Wilmington/CMS Work/Summer (0516-0728)/Main Functions')
addpath('C:\Users\nsc4750\OneDrive - UNC-Wilmington\CMS Work\Summer (0516-0728)\Main Functions') 
%addpath('../../Summer (0516-0728)/Main Function')

%% Calculate
% - Do all of the energy weighted means per frequency band and save
%   everything to varibales to be plotted later on.


Bnames = fields(Buoy.Obs); % to list buoy names


for MO = 1:5 % for observed and 4 models
    
    for i = [1:4 6:8] % for the 9 buoys (7 pairs)
        
        if MO == 1
            eval(sprintf('STODO1 = Buoy.Obs.%s;',Bnames{i})) % Observed
            eval(sprintf('STODO2 = Buoy.Obs.%s;',Bnames{i+1}))
        elseif MO == 2
            eval(sprintf('STODO1 = Buoy.SWAN050.%s;',Bnames{i})) % SWAN050
            eval(sprintf('STODO2 = Buoy.SWAN050.%s;',Bnames{i+1}))
        elseif MO == 3
            eval(sprintf('STODO1 = Buoy.SWAN050sm.%s;',Bnames{i})) % SWAN050sm
            eval(sprintf('STODO2 = Buoy.SWAN050sm.%s;',Bnames{i+1}))
        elseif MO == 4
            eval(sprintf('STODO1 = Buoy.COUP050.%s;',Bnames{i})) % COUP050
            eval(sprintf('STODO2 = Buoy.COUP050.%s;',Bnames{i+1}))
        else
            eval(sprintf('STODO1 = Buoy.SWAN075.%s;',Bnames{i})) % SWAN075
            eval(sprintf('STODO2 = Buoy.SWAN075.%s;',Bnames{i+1}))
        end
        
        
        
        for j = 1:length(Buoy.time) % for each individual hour
            fdivfp1 = Buoy.freq./STODO1.pfreq(j);
            fdivfp2 = Buoy.freq./STODO2.pfreq(j);
            fdivfp = (fdivfp1 + fdivfp2)./2;
            
            ind = find(fdivfp>0.5 & fdivfp<=0.8);
            [TED,~,~,~] = NC_ObsDissDF(STODO1.See(ind,j),STODO2.See(ind,j),...
                STODO1.Dir(ind,j),STODO2.Dir(ind,j),Buoy.freq(ind),...
                STODO1.utm,STODO2.utm,STODO1.depth(j),STODO2.depth(j),0.01);
            avgSee = (STODO1.See(ind,j) + STODO2.See(ind,j))./2;
            Am0 = sum(avgSee).*0.01;
            m1TED = sum(avgSee.*TED).*0.01;
            EWM_TED1(j) = m1TED./Am0;
            if EWM_TED1(j) > 40 || EWM_TED1(j) < -40 % remove outliers
                EWM_TED1(j) = NaN;
            end
            diffD = STODO2.Dir(ind,j) - STODO1.Dir(ind,j); 
            Dm1 = sum(diffD.*avgSee).*0.01;
            D_EWM1(j) = Dm1./Am0;
            if D_EWM1(j) < -50 % remove outliers
                D_EWM1(j) = NaN;
            end
            diffDS = STODO2.Dspread(ind,j) - STODO1.Dspread(ind,j); 
            DSm1 = sum(diffDS.*avgSee).*0.01;
            DS_EWM1(j) = DSm1./Am0;
            
            
            ind = find(fdivfp>0.8 & fdivfp<=1.2);
            [TED,~,~,~] = NC_ObsDissDF(STODO1.See(ind,j),STODO2.See(ind,j),...
                STODO1.Dir(ind,j),STODO2.Dir(ind,j),Buoy.freq(ind),...
                STODO1.utm,STODO2.utm,STODO1.depth(j),STODO2.depth(j),0.01);
            avgSee = (STODO1.See(ind,j) + STODO2.See(ind,j))./2;
            Am0 = sum(avgSee).*0.01;
            m1TED = sum(avgSee.*TED).*0.01;
            EWM_TED2(j) = m1TED./Am0;
            if EWM_TED2(j) > 40 || EWM_TED2(j) < -40 % remove outliers
                EWM_TED2(j) = NaN;
            end
            diffD = STODO2.Dir(ind,j) - STODO1.Dir(ind,j); 
            Dm1 = sum(diffD.*avgSee).*0.01;
            D_EWM2(j) = Dm1./Am0;
            if D_EWM2(j) < -50 % remove outliers
                D_EWM2(j) = NaN;
            end
            diffDS = STODO2.Dspread(ind,j) - STODO1.Dspread(ind,j); 
            DSm1 = sum(diffDS.*avgSee).*0.01;
            DS_EWM2(j) = DSm1./Am0;

            
            ind = find(fdivfp>1.2 & fdivfp<=2.0);
            [TED,~,~,~] = NC_ObsDissDF(STODO1.See(ind,j),STODO2.See(ind,j),...
                STODO1.Dir(ind,j),STODO2.Dir(ind,j),Buoy.freq(ind),...
                STODO1.utm,STODO2.utm,STODO1.depth(j),STODO2.depth(j),0.01);
            avgSee = (STODO1.See(ind,j) + STODO2.See(ind,j))./2;
            Am0 = sum(avgSee).*0.01;
            m1TED = sum(avgSee.*TED).*0.01;
            EWM_TED3(j) = m1TED./Am0;
            if EWM_TED3(j) > 40 || EWM_TED3(j) < -40 % remove outliers
                EWM_TED3(j) = NaN;
            end
            diffD = STODO2.Dir(ind,j) - STODO1.Dir(ind,j); 
            Dm1 = sum(diffD.*avgSee).*0.01;
            D_EWM3(j) = Dm1./Am0;
            if D_EWM3(j) < -50 % remove outliers
                D_EWM3(j) = NaN;
            end
            diffDS = STODO2.Dspread(ind,j) - STODO1.Dspread(ind,j); 
            DSm1 = sum(diffDS.*avgSee).*0.01;
            DS_EWM3(j) = DSm1./Am0;

            
            ind = find(fdivfp>2.0 & fdivfp<=3.0);
            [TED,~,~,~] = NC_ObsDissDF(STODO1.See(ind,j),STODO2.See(ind,j),...
                STODO1.Dir(ind,j),STODO2.Dir(ind,j),Buoy.freq(ind),...
                STODO1.utm,STODO2.utm,STODO1.depth(j),STODO2.depth(j),0.01);
            avgSee = (STODO1.See(ind,j) + STODO2.See(ind,j))./2;
            Am0 = sum(avgSee).*0.01;
            m1TED = sum(avgSee.*TED).*0.01;
            EWM_TED4(j) = m1TED./Am0;
            if EWM_TED4(j) > 5 || EWM_TED4(j) < -5 % remove outliers
                EWM_TED4(j) = NaN;
            end
            diffD = STODO2.Dir(ind,j) - STODO1.Dir(ind,j); 
            Dm1 = sum(diffD.*avgSee).*0.01;
            D_EWM4(j) = Dm1./Am0;
            if D_EWM4(j) < -100 % remove outliers
                D_EWM4(j) = NaN;
            end
            diffDS = STODO2.Dspread(ind,j) - STODO1.Dspread(ind,j); 
            DSm1 = sum(diffDS.*avgSee).*0.01;
            DS_EWM4(j) = DSm1./Am0;
            
        end
        
        if MO == 1
            eval(sprintf('Buoy.Obs.%s.EWM_TED1 = EWM_TED1;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWM_TED2 = EWM_TED2;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWM_TED3 = EWM_TED3;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWM_TED4 = EWM_TED4;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMdiffDir1 = D_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMdiffDir2 = D_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMdiffDir3 = D_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMdiffDir4 = D_EWM4;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMdiffDspread1 = DS_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMdiffDspread2 = DS_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMdiffDspread3 = DS_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMdiffDspread4 = DS_EWM4;',Bnames{i}))
        elseif MO == 2
            eval(sprintf('Buoy.SWAN050.%s.EWM_TED1 = EWM_TED1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWM_TED2 = EWM_TED2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWM_TED3 = EWM_TED3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWM_TED4 = EWM_TED4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMdiffDir1 = D_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMdiffDir2 = D_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMdiffDir3 = D_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMdiffDir4 = D_EWM4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMdiffDspread1 = DS_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMdiffDspread2 = DS_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMdiffDspread3 = DS_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMdiffDspread4 = DS_EWM4;',Bnames{i}))
        elseif MO == 3
            eval(sprintf('Buoy.SWAN050sm.%s.EWM_TED1 = EWM_TED1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWM_TED2 = EWM_TED2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWM_TED3 = EWM_TED3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWM_TED4 = EWM_TED4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMdiffDir1 = D_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMdiffDir2 = D_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMdiffDir3 = D_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMdiffDir4 = D_EWM4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMdiffDspread1 = DS_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMdiffDspread2 = DS_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMdiffDspread3 = DS_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMdiffDspread4 = DS_EWM4;',Bnames{i}))
        elseif MO ==4
            eval(sprintf('Buoy.COUP050.%s.EWM_TED1 = EWM_TED1;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWM_TED2 = EWM_TED2;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWM_TED3 = EWM_TED3;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWM_TED4 = EWM_TED4;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMdiffDir1 = D_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMdiffDir2 = D_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMdiffDir3 = D_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMdiffDir4 = D_EWM4;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMdiffDspread1 = DS_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMdiffDspread2 = DS_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMdiffDspread3 = DS_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMdiffDspread4 = DS_EWM4;',Bnames{i}))
        else
            eval(sprintf('Buoy.SWAN075.%s.EWM_TED1 = EWM_TED1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWM_TED2 = EWM_TED2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWM_TED3 = EWM_TED3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWM_TED4 = EWM_TED4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMdiffDir1 = D_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMdiffDir2 = D_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMdiffDir3 = D_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMdiffDir4 = D_EWM4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMdiffDspread1 = DS_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMdiffDspread2 = DS_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMdiffDspread3 = DS_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMdiffDspread4 = DS_EWM4;',Bnames{i}))
        end
        
    end
end


%% Plot Formatting

lightgrey = [0.85 0.85 0.85]; %the color of the background of subplots
pformat = {'xb','xm','xr','xc','xk','ob','om','or','oc'}; % The color and shape of markers
transformat = [1:-1/11:0.27]; %the transparency of the markers
MSize = 35; %MarkerSize

one_to_one = -30:360; % to plot a one to one line through plots

BTlegendnames = {'B01-B03','B03-B05','B05-B10','B10-B13','X01-X03','X03-X04','X04-X05'};
%BTlegendnamesS = {'B01-B03','B03-B05'};


%% Plot TED Comparison

figure(4);clf; % Energy Figure
% set(gcf,'position',[600,89,1282,863])
set(gcf,'position',[600,50,1000,940])

    % SWAN050 Comparison:
    
SWAN_TEDO{1} = []; SWAN_TEDM{1} = [];
SWAN_TEDO{2} = []; SWAN_TEDM{2} = [];
SWAN_TEDO{3} = []; SWAN_TEDM{3} = [];
SWAN_TEDO{4} = []; SWAN_TEDM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OTED1 = Buoy.Obs.%s.EWM_TED1;',Bnames{i}))
    eval(sprintf('MTED1 = Buoy.SWAN050.%s.EWM_TED1;',Bnames{i}))
    eval(sprintf('OTED2 = Buoy.Obs.%s.EWM_TED2;',Bnames{i}))
    eval(sprintf('MTED2 = Buoy.SWAN050.%s.EWM_TED2;',Bnames{i}))
    eval(sprintf('OTED3 = Buoy.Obs.%s.EWM_TED3;',Bnames{i}))
    eval(sprintf('MTED3 = Buoy.SWAN050.%s.EWM_TED3;',Bnames{i}))
    eval(sprintf('OTED4 = Buoy.Obs.%s.EWM_TED4;',Bnames{i}))
    eval(sprintf('MTED4 = Buoy.SWAN050.%s.EWM_TED4;',Bnames{i}))
    
    subplot(441)
    scatter(OTED1,MTED1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(442)
    scatter(OTED2,MTED2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(443)
    scatter(OTED3,MTED3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(444)
    scatter(OTED4,MTED4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN_TEDO{1} = cat(2,SWAN_TEDO{1},OTED1);
    SWAN_TEDO{2} = cat(2,SWAN_TEDO{2},OTED2);
    SWAN_TEDO{3} = cat(2,SWAN_TEDO{3},OTED3);
    SWAN_TEDO{4} = cat(2,SWAN_TEDO{4},OTED4);
    
    SWAN_TEDM{1} = cat(2,SWAN_TEDM{1},MTED1);
    SWAN_TEDM{2} = cat(2,SWAN_TEDM{2},MTED2);
    SWAN_TEDM{3} = cat(2,SWAN_TEDM{3},MTED3);
    SWAN_TEDM{4} = cat(2,SWAN_TEDM{4},MTED4);
end

subplot(441)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.5>f/f_p\geq0.8)')
xlabel('Obs TED (W/m)')
ylabel('SWAN050 TED (W/m)')
legend(BTlegendnames,'NumColumns',2,'location',[0.418762849685225,0.003724780584812,0.185999997496605,0.065425530139436])
[~,Sk{1},maxcor{1},rmse{1}] = mod_error(SWAN_TEDM{1},SWAN_TEDO{1});
STD{1} = nanstd(SWAN_TEDM{1});
OSTD{1} = nanstd(SWAN_TEDO{1});
str = {sprintf('RMSE=%.4f',round(rmse{1}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{1}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{1}.CCall(101),4))};
dim = [0.132897035881435,0.871857202731689,0.048190327613105,0.046384704519119];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(442)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.8>f/f_p\geq1.2)')
xlabel('Obs TED (W/m)')
ylabel('SWAN050 TED (W/m)')
[~,Sk{2},maxcor{2},rmse{2}] = mod_error(SWAN_TEDM{2},SWAN_TEDO{2});
STD{2} = nanstd(SWAN_TEDM{2});
OSTD{2} = nanstd(SWAN_TEDO{2});
str = {sprintf('RMSE=%.4f',round(rmse{2}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{2}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{2}.CCall(101),4))};
dim = [0.337188560639638,0.853992217708371,0.067445363700615,0.063882062344352];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(443)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (1.2>f/f_p\geq2.0)')
xlabel('Obs TED (W/m)')
ylabel('SWAN050 TED (W/m)')
[~,Sk{3},maxcor{3},rmse{3}] = mod_error(SWAN_TEDM{3},SWAN_TEDO{3});
STD{3} = nanstd(SWAN_TEDM{3});
OSTD{3} = nanstd(SWAN_TEDO{3});
str = {sprintf('RMSE=%.4f',round(rmse{3}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{3}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{3}.CCall(101),4))};
dim = [0.548205928237131,0.855301859718654,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(444)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (2.0>f/f_p\geq3.0)')
xlabel('Obs TED (W/m)')
ylabel('SWAN050 TED (W/m)')
[~,Sk{4},maxcor{4},rmse{4}] = mod_error(SWAN_TEDM{4},SWAN_TEDO{4});
STD{4} = nanstd(SWAN_TEDM{4});
OSTD{4} = nanstd(SWAN_TEDO{4});
str = {sprintf('RMSE=%.4f',round(rmse{4}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{4}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{4}.CCall(101),4))};
dim = [0.750390015600627,0.856460608270219,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN050sm Comparison:
    
SWANsm_TEDO{1} = []; SWANsm_TEDM{1} = [];
SWANsm_TEDO{2} = []; SWANsm_TEDM{2} = [];
SWANsm_TEDO{3} = []; SWANsm_TEDM{3} = [];
SWANsm_TEDO{4} = []; SWANsm_TEDM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OTED1 = Buoy.Obs.%s.EWM_TED1;',Bnames{i}))
    eval(sprintf('MTED1 = Buoy.SWAN050sm.%s.EWM_TED1;',Bnames{i}))
    eval(sprintf('OTED2 = Buoy.Obs.%s.EWM_TED2;',Bnames{i}))
    eval(sprintf('MTED2 = Buoy.SWAN050sm.%s.EWM_TED2;',Bnames{i}))
    eval(sprintf('OTED3 = Buoy.Obs.%s.EWM_TED3;',Bnames{i}))
    eval(sprintf('MTED3 = Buoy.SWAN050sm.%s.EWM_TED3;',Bnames{i}))
    eval(sprintf('OTED4 = Buoy.Obs.%s.EWM_TED4;',Bnames{i}))
    eval(sprintf('MTED4 = Buoy.SWAN050sm.%s.EWM_TED4;',Bnames{i}))
    
    subplot(445)
    scatter(OTED1,MTED1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(446)
    scatter(OTED2,MTED2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(447)
    scatter(OTED3,MTED3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(448)
    scatter(OTED4,MTED4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWANsm_TEDO{1} = cat(2,SWANsm_TEDO{1},OTED1);
    SWANsm_TEDO{2} = cat(2,SWANsm_TEDO{2},OTED2);
    SWANsm_TEDO{3} = cat(2,SWANsm_TEDO{3},OTED3);
    SWANsm_TEDO{4} = cat(2,SWANsm_TEDO{4},OTED4);
    
    SWANsm_TEDM{1} = cat(2,SWANsm_TEDM{1},MTED1);
    SWANsm_TEDM{2} = cat(2,SWANsm_TEDM{2},MTED2);
    SWANsm_TEDM{3} = cat(2,SWANsm_TEDM{3},MTED3);
    SWANsm_TEDM{4} = cat(2,SWANsm_TEDM{4},MTED4);
end

subplot(445)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.5>f/f_p\geq0.8)')
xlabel('Obs TED (W/m)')
ylabel('SWAN050 TED (W/m)')
%legend(BTlegendnames,'NumColumns',2,'location','north')
[~,Sk{5},maxcor{5},rmse{5}] = mod_error(SWANsm_TEDM{1},SWANsm_TEDO{1});
STD{5} = nanstd(SWANsm_TEDM{1});
str = {sprintf('RMSE=%.4f',round(rmse{5}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{5}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{5}.CCall(101),4))};
dim = [0.1337,0.6378,0.0698,0.0603];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(446)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.8>f/f_p\geq1.2)')
xlabel('Obs TED (W/m)')
ylabel('SWAN050 TED (W/m)')
[~,Sk{6},maxcor{6},rmse{6}] = mod_error(SWANsm_TEDM{2},SWANsm_TEDO{2});
STD{6} = nanstd(SWANsm_TEDM{2});
str = {sprintf('RMSE=%.4f',round(rmse{6}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{6}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{6}.CCall(101),4))};
dim = [0.340873634945401,0.639774629127692,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(447)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (1.2>f/f_p\geq2.0)')
xlabel('Obs TED (W/m)')
ylabel('SWAN050 TED (W/m)')
[~,Sk{7},maxcor{7},rmse{7}] = mod_error(SWANsm_TEDM{3},SWANsm_TEDO{3});
STD{7} = nanstd(SWANsm_TEDM{3});
str = {sprintf('RMSE=%.4f',round(rmse{7}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{7}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{7}.CCall(101),4))};
dim = [0.548845553822153,0.632363848495593,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(448)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (2.0>f/f_p\geq3.0)')
xlabel('Obs TED (W/m)')
ylabel('SWAN050 TED (W/m)')
[~,Sk{8},maxcor{8},rmse{8}] = mod_error(SWANsm_TEDM{4},SWANsm_TEDO{4});
STD{8} = nanstd(SWANsm_TEDM{4});
str = {sprintf('RMSE=%.4f',round(rmse{8}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{8}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{8}.CCall(101),4))};
dim = [0.753213728549143,0.635840094150286,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % COUP050 Comparison:
    
COUP_TEDO{1} = []; COUP_TEDM{1} = [];
COUP_TEDO{2} = []; COUP_TEDM{2} = [];
COUP_TEDO{3} = []; COUP_TEDM{3} = [];
COUP_TEDO{4} = []; COUP_TEDM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OTED1 = Buoy.Obs.%s.EWM_TED1;',Bnames{i}))
    eval(sprintf('MTED1 = Buoy.COUP050.%s.EWM_TED1;',Bnames{i}))
    eval(sprintf('OTED2 = Buoy.Obs.%s.EWM_TED2;',Bnames{i}))
    eval(sprintf('MTED2 = Buoy.COUP050.%s.EWM_TED2;',Bnames{i}))
    eval(sprintf('OTED3 = Buoy.Obs.%s.EWM_TED3;',Bnames{i}))
    eval(sprintf('MTED3 = Buoy.COUP050.%s.EWM_TED3;',Bnames{i}))
    eval(sprintf('OTED4 = Buoy.Obs.%s.EWM_TED4;',Bnames{i}))
    eval(sprintf('MTED4 = Buoy.COUP050.%s.EWM_TED4;',Bnames{i}))
    
    subplot(449)
    scatter(OTED1,MTED1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,10)
    scatter(OTED2,MTED2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,11)
    scatter(OTED3,MTED3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,12)
    scatter(OTED4,MTED4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    COUP_TEDO{1} = cat(2,COUP_TEDO{1},OTED1);
    COUP_TEDO{2} = cat(2,COUP_TEDO{2},OTED2);
    COUP_TEDO{3} = cat(2,COUP_TEDO{3},OTED3);
    COUP_TEDO{4} = cat(2,COUP_TEDO{4},OTED4);
    
    COUP_TEDM{1} = cat(2,COUP_TEDM{1},MTED1);
    COUP_TEDM{2} = cat(2,COUP_TEDM{2},MTED2);
    COUP_TEDM{3} = cat(2,COUP_TEDM{3},MTED3);
    COUP_TEDM{4} = cat(2,COUP_TEDM{4},MTED4);
end

subplot(449)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.5>f/f_p\geq0.8)')
xlabel('Obs TED (W/m)')
ylabel('COUP050 TED (W/m)')
%legend(BTlegendnames,'NumColumns',2,'location','north')
[~,Sk{9},maxcor{9},rmse{9}] = mod_error(COUP_TEDM{1},COUP_TEDO{1});
STD{9} = nanstd(COUP_TEDM{1});
str = {sprintf('RMSE=%.4f',round(rmse{9}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{9}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{9}.CCall(101),4))};
dim = [0.132868954758192,0.419632406963086,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,10)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.8>f/f_p\geq1.2)')
xlabel('Obs TED (W/m)')
ylabel('COUP050 TED (W/m)')
[~,Sk{10},maxcor{10},rmse{10}] = mod_error(COUP_TEDM{2},COUP_TEDO{2});
STD{10} = nanstd(COUP_TEDM{2});
str = {sprintf('RMSE=%.4f',round(rmse{10}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{10}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{10}.CCall(101),4))};
dim = [0.338237129485182,0.417995366456196,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,11)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (1.2>f/f_p\geq2.0)')
xlabel('Obs TED (W/m)')
ylabel('COUP050 TED (W/m)')
[~,Sk{11},maxcor{11},rmse{11}] = mod_error(COUP_TEDM{3},COUP_TEDO{3});
STD{11} = nanstd(COUP_TEDM{3});
str = {sprintf('RMSE=%.4f',round(rmse{11}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{11}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{11}.CCall(101),4))};
dim = [0.544945397815915,0.417995366456196,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,12)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (2.0>f/f_p\geq3.0)')
xlabel('Obs TED (W/m)')
ylabel('COUP050 TED (W/m)')
[~,Sk{12},maxcor{12},rmse{12}] = mod_error(COUP_TEDM{4},COUP_TEDO{4});
STD{12} = nanstd(COUP_TEDM{4});
str = {sprintf('RMSE=%.4f',round(rmse{12}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{12}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{12}.CCall(101),4))};
dim = [0.750468018720751,0.417647741890727,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN075 Comparison:
    
SWAN075_TEDO{1} = []; SWAN075_TEDM{1} = [];
SWAN075_TEDO{2} = []; SWAN075_TEDM{2} = [];
SWAN075_TEDO{3} = []; SWAN075_TEDM{3} = [];
SWAN075_TEDO{4} = []; SWAN075_TEDM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OTED1 = Buoy.Obs.%s.EWM_TED1;',Bnames{i}))
    eval(sprintf('MTED1 = Buoy.SWAN075.%s.EWM_TED1;',Bnames{i}))
    eval(sprintf('OTED2 = Buoy.Obs.%s.EWM_TED2;',Bnames{i}))
    eval(sprintf('MTED2 = Buoy.SWAN075.%s.EWM_TED2;',Bnames{i}))
    eval(sprintf('OTED3 = Buoy.Obs.%s.EWM_TED3;',Bnames{i}))
    eval(sprintf('MTED3 = Buoy.SWAN075.%s.EWM_TED3;',Bnames{i}))
    eval(sprintf('OTED4 = Buoy.Obs.%s.EWM_TED4;',Bnames{i}))
    eval(sprintf('MTED4 = Buoy.SWAN075.%s.EWM_TED4;',Bnames{i}))
    
    subplot(4,4,13)
    scatter(OTED1,MTED1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,14)
    scatter(OTED2,MTED2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,15)
    scatter(OTED3,MTED3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,16)
    scatter(OTED4,MTED4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN075_TEDO{1} = cat(2,SWAN075_TEDO{1},OTED1);
    SWAN075_TEDO{2} = cat(2,SWAN075_TEDO{2},OTED2);
    SWAN075_TEDO{3} = cat(2,SWAN075_TEDO{3},OTED3);
    SWAN075_TEDO{4} = cat(2,SWAN075_TEDO{4},OTED4);
    
    SWAN075_TEDM{1} = cat(2,SWAN075_TEDM{1},MTED1);
    SWAN075_TEDM{2} = cat(2,SWAN075_TEDM{2},MTED2);
    SWAN075_TEDM{3} = cat(2,SWAN075_TEDM{3},MTED3);
    SWAN075_TEDM{4} = cat(2,SWAN075_TEDM{4},MTED4);
end

subplot(4,4,13)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.5>f/f_p\geq0.8)')
xlabel('Obs TED (W/m)')
ylabel('SWAN075 TED (W/m)')
%legend(BTlegendnames,'NumColumns',2,'location','north')
[~,Sk{13},maxcor{13},rmse{13}] = mod_error(SWAN075_TEDM{1},SWAN075_TEDO{1});
STD{13} = nanstd(SWAN075_TEDM{1});
str = {sprintf('RMSE=%.4f',round(rmse{13}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{13}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{13}.CCall(101),4))};
dim = [0.1334,0.2028,0.0698,0.0603];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,14)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.8>f/f_p\geq1.2)')
xlabel('Obs TED (W/m)')
ylabel('SWAN075 TED (W/m)')
[~,Sk{14},maxcor{14},rmse{14}] = mod_error(SWAN075_TEDM{2},SWAN075_TEDO{2});
STD{14} = nanstd(SWAN075_TEDM{2});
str = {sprintf('RMSE=%.4f',round(rmse{14}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{14}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{14}.CCall(101),4))};
dim = [0.339797191887678,0.197833141658977,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,15)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (1.2>f/f_p\geq2.0)')
xlabel('Obs TED (W/m)')
ylabel('SWAN075 TED (W/m)')
[~,Sk{15},maxcor{15},rmse{15}] = mod_error(SWAN075_TEDM{3},SWAN075_TEDO{3});
STD{15} = nanstd(SWAN075_TEDM{3});
str = {sprintf('RMSE=%.4f',round(rmse{15}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{15}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{15}.CCall(101),4))};
dim = [0.547285491419659,0.196674393107413,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,16)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (2.0>f/f_p\geq3.0)')
xlabel('Obs TED (W/m)')
ylabel('SWAN075 TED (W/m)')
[~,Sk{16},maxcor{16},rmse{16}] = mod_error(SWAN075_TEDM{4},SWAN075_TEDO{4});
STD{16} = nanstd(SWAN075_TEDM{4});
str = {sprintf('RMSE=%.4f',round(rmse{16}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{16}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{16}.CCall(101),4))};
dim = [0.751248049921999,0.197485517093508,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')



% Rename TED Statistics Variables
TED_rmse = rmse;
TED_Sk = Sk;
TED_maxcor = maxcor;
TED_STD = STD;
TED_OSTD = OSTD;
clear rmse Sk maxcor STD OSTD


%% Plot Difference in Direction b/t Buoys

figure(5);clf; % Energy Figure
set(gcf,'position',[750,50,1000,940])

    % SWAN050 Comparison: 
SWAN_diffDO{1} = []; SWAN_diffDM{1} = [];
SWAN_diffDO{2} = []; SWAN_diffDM{2} = [];
SWAN_diffDO{3} = []; SWAN_diffDM{3} = [];
SWAN_diffDO{4} = []; SWAN_diffDM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OdiffDir1 = Buoy.Obs.%s.EWMdiffDir1;',Bnames{i}))
    eval(sprintf('MdiffDir1 = Buoy.SWAN050.%s.EWMdiffDir1;',Bnames{i}))
    eval(sprintf('OdiffDir2 = Buoy.Obs.%s.EWMdiffDir2;',Bnames{i}))
    eval(sprintf('MdiffDir2 = Buoy.SWAN050.%s.EWMdiffDir2;',Bnames{i}))
    eval(sprintf('OdiffDir3 = Buoy.Obs.%s.EWMdiffDir3;',Bnames{i}))
    eval(sprintf('MdiffDir3 = Buoy.SWAN050.%s.EWMdiffDir3;',Bnames{i}))
    eval(sprintf('OdiffDir4 = Buoy.Obs.%s.EWMdiffDir4;',Bnames{i}))
    eval(sprintf('MdiffDir4 = Buoy.SWAN050.%s.EWMdiffDir4;',Bnames{i}))
    
    subplot(441)
    scatter(OdiffDir1,MdiffDir1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(442)
    scatter(OdiffDir2,MdiffDir2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(443)
    scatter(OdiffDir3,MdiffDir3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(444)
    scatter(OdiffDir4,MdiffDir4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN_diffDO{1} = cat(2,SWAN_diffDO{1},OdiffDir1);
    SWAN_diffDO{2} = cat(2,SWAN_diffDO{2},OdiffDir2);
    SWAN_diffDO{3} = cat(2,SWAN_diffDO{3},OdiffDir3);
    SWAN_diffDO{4} = cat(2,SWAN_diffDO{4},OdiffDir4);
    
    SWAN_diffDM{1} = cat(2,SWAN_diffDM{1},MdiffDir1);
    SWAN_diffDM{2} = cat(2,SWAN_diffDM{2},MdiffDir2);
    SWAN_diffDM{3} = cat(2,SWAN_diffDM{3},MdiffDir3);
    SWAN_diffDM{4} = cat(2,SWAN_diffDM{4},MdiffDir4);
end

subplot(441)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.5>f/f_p\geq0.8)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN050 diffDir (^o)')
legend(BTlegendnames,'NumColumns',2,'location',[0.418762849685225,0.003724780584812,0.185999997496605,0.065425530139436])
[~,Sk{1},maxcor{1},rmse{1}] = mod_error(SWAN_diffDM{1},SWAN_diffDO{1});
STD{1} = nanstd(SWAN_diffDM{1});
OSTD{1} = nanstd(SWAN_diffDO{1});
str = {sprintf('RMSE=%.4f',round(rmse{1}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{1}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{1}.CCall(101),4))};
dim = [0.1407,0.8522,0.0674,0.0639];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(442)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.8>f/f_p\geq1.2)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN050 diffDir (^o)')
[~,Sk{2},maxcor{2},rmse{2}] = mod_error(SWAN_diffDM{2},SWAN_diffDO{2});
STD{2} = nanstd(SWAN_diffDM{2});
OSTD{2} = nanstd(SWAN_diffDO{2});
str = {sprintf('RMSE=%.4f',round(rmse{2}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{2}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{2}.CCall(101),4))};
dim = [0.339528654243382,0.852833469156808,0.067445363700615,0.063882062344352];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(443)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (1.2>f/f_p\geq2.0)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN050 diffDir (^o)')
[~,Sk{3},maxcor{3},rmse{3}] = mod_error(SWAN_diffDM{3},SWAN_diffDO{3});
STD{3} = nanstd(SWAN_diffDM{3});
OSTD{3} = nanstd(SWAN_diffDO{3});
str = {sprintf('RMSE=%.4f',round(rmse{3}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{3}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{3}.CCall(101),4))};
dim = [0.545456891372867,0.855150966259937,0.067445363700615,0.063882062344352];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(444)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (2.0>f/f_p\geq3.0)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN050 diffDir (^o)')
[~,Sk{4},maxcor{4},rmse{4}] = mod_error(SWAN_diffDM{4},SWAN_diffDO{4});
STD{4} = nanstd(SWAN_diffDM{4});
OSTD{4} = nanstd(SWAN_diffDO{4});
str = {sprintf('RMSE=%.4f',round(rmse{4}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{4}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{4}.CCall(101),4))};
dim = [0.751950078003123,0.852984362615525,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN050sm Comparison: 
SWANsm_diffDO{1} = []; SWANsm_diffDM{1} = [];
SWANsm_diffDO{2} = []; SWANsm_diffDM{2} = [];
SWANsm_diffDO{3} = []; SWANsm_diffDM{3} = [];
SWANsm_diffDO{4} = []; SWANsm_diffDM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OdiffDir1 = Buoy.Obs.%s.EWMdiffDir1;',Bnames{i}))
    eval(sprintf('MdiffDir1 = Buoy.SWAN050sm.%s.EWMdiffDir1;',Bnames{i}))
    eval(sprintf('OdiffDir2 = Buoy.Obs.%s.EWMdiffDir2;',Bnames{i}))
    eval(sprintf('MdiffDir2 = Buoy.SWAN050sm.%s.EWMdiffDir2;',Bnames{i}))
    eval(sprintf('OdiffDir3 = Buoy.Obs.%s.EWMdiffDir3;',Bnames{i}))
    eval(sprintf('MdiffDir3 = Buoy.SWAN050sm.%s.EWMdiffDir3;',Bnames{i}))
    eval(sprintf('OdiffDir4 = Buoy.Obs.%s.EWMdiffDir4;',Bnames{i}))
    eval(sprintf('MdiffDir4 = Buoy.SWAN050sm.%s.EWMdiffDir4;',Bnames{i}))
    
    subplot(445)
    scatter(OdiffDir1,MdiffDir1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(446)
    scatter(OdiffDir2,MdiffDir2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(447)
    scatter(OdiffDir3,MdiffDir3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(448)
    scatter(OdiffDir4,MdiffDir4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWANsm_diffDO{1} = cat(2,SWANsm_diffDO{1},OdiffDir1);
    SWANsm_diffDO{2} = cat(2,SWANsm_diffDO{2},OdiffDir2);
    SWANsm_diffDO{3} = cat(2,SWANsm_diffDO{3},OdiffDir3);
    SWANsm_diffDO{4} = cat(2,SWANsm_diffDO{4},OdiffDir4);
    
    SWANsm_diffDM{1} = cat(2,SWANsm_diffDM{1},MdiffDir1);
    SWANsm_diffDM{2} = cat(2,SWANsm_diffDM{2},MdiffDir2);
    SWANsm_diffDM{3} = cat(2,SWANsm_diffDM{3},MdiffDir3);
    SWANsm_diffDM{4} = cat(2,SWANsm_diffDM{4},MdiffDir4);
end

subplot(445)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.5>f/f_p\geq0.8)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN050sm diffDir (^o)')
%legend(BTlegendnames,'NumColumns',2,'location','north')
[~,Sk{5},maxcor{5},rmse{5}] = mod_error(SWANsm_diffDM{1},SWANsm_diffDO{1});
STD{5} = nanstd(SWANsm_diffDM{1});
str = {sprintf('RMSE=%.4f',round(rmse{5}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{5}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{5}.CCall(101),4))};
dim = [0.1375,0.6392,0.0698,0.0603];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(446)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.8>f/f_p\geq1.2)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN050sm diffDir (^o)')
[~,Sk{6},maxcor{6},rmse{6}] = mod_error(SWANsm_diffDM{2},SWANsm_diffDO{2});
STD{6} = nanstd(SWANsm_diffDM{2});
str = {sprintf('RMSE=%.4f',round(rmse{6}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{6}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{6}.CCall(101),4))};
dim = [0.338533541341657,0.636298383472999,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(447)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (1.2>f/f_p\geq2.0)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN050sm diffDir (^o)')
[~,Sk{7},maxcor{7},rmse{7}] = mod_error(SWANsm_diffDM{3},SWANsm_diffDO{3});
STD{7} = nanstd(SWANsm_diffDM{3});
str = {sprintf('RMSE=%.4f',round(rmse{7}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{7}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{7}.CCall(101),4))};
dim = [0.54650546021841,0.639316339804978,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(448)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (2.0>f/f_p\geq3.0)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN050sm diffDir (^o)')
[~,Sk{8},maxcor{8},rmse{8}] = mod_error(SWANsm_diffDM{4},SWANsm_diffDO{4});
STD{8} = nanstd(SWANsm_diffDM{4});
str = {sprintf('RMSE=%.4f',round(rmse{8}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{8}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{8}.CCall(101),4))};
dim = [0.752589703588144,0.638157591253414,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % COUP050 Comparison: 
COUP050_diffDO{1} = []; COUP050_diffDM{1} = [];
COUP050_diffDO{2} = []; COUP050_diffDM{2} = [];
COUP050_diffDO{3} = []; COUP050_diffDM{3} = [];
COUP050_diffDO{4} = []; COUP050_diffDM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OdiffDir1 = Buoy.Obs.%s.EWMdiffDir1;',Bnames{i}))
    eval(sprintf('MdiffDir1 = Buoy.COUP050.%s.EWMdiffDir1;',Bnames{i}))
    eval(sprintf('OdiffDir2 = Buoy.Obs.%s.EWMdiffDir2;',Bnames{i}))
    eval(sprintf('MdiffDir2 = Buoy.COUP050.%s.EWMdiffDir2;',Bnames{i}))
    eval(sprintf('OdiffDir3 = Buoy.Obs.%s.EWMdiffDir3;',Bnames{i}))
    eval(sprintf('MdiffDir3 = Buoy.COUP050.%s.EWMdiffDir3;',Bnames{i}))
    eval(sprintf('OdiffDir4 = Buoy.Obs.%s.EWMdiffDir4;',Bnames{i}))
    eval(sprintf('MdiffDir4 = Buoy.COUP050.%s.EWMdiffDir4;',Bnames{i}))
    
    subplot(449)
    scatter(OdiffDir1,MdiffDir1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,10)
    scatter(OdiffDir2,MdiffDir2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,11)
    scatter(OdiffDir3,MdiffDir3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,12)
    scatter(OdiffDir4,MdiffDir4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    COUP050_diffDO{1} = cat(2,COUP050_diffDO{1},OdiffDir1);
    COUP050_diffDO{2} = cat(2,COUP050_diffDO{2},OdiffDir2);
    COUP050_diffDO{3} = cat(2,COUP050_diffDO{3},OdiffDir3);
    COUP050_diffDO{4} = cat(2,COUP050_diffDO{4},OdiffDir4);
    
    COUP050_diffDM{1} = cat(2,COUP050_diffDM{1},MdiffDir1);
    COUP050_diffDM{2} = cat(2,COUP050_diffDM{2},MdiffDir2);
    COUP050_diffDM{3} = cat(2,COUP050_diffDM{3},MdiffDir3);
    COUP050_diffDM{4} = cat(2,COUP050_diffDM{4},MdiffDir4);
end

subplot(449)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.5>f/f_p\geq0.8)')
xlabel('Obs diffDir (^o)')
ylabel('COUP050 diffDir (^o)')
%legend(BTlegendnames,'NumColumns',2,'location','north')
[~,Sk{9},maxcor{9},rmse{9}] = mod_error(COUP050_diffDM{1},COUP050_diffDO{1});
STD{9} = nanstd(COUP050_diffDM{1});
str = {sprintf('RMSE=%.4f',round(rmse{9}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{9}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{9}.CCall(101),4))};
dim = [0.139301092043682,0.41598234902566,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,10)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.8>f/f_p\geq1.2)')
xlabel('Obs diffDir (^o)')
ylabel('COUP050 diffDir (^o)')
[~,Sk{10},maxcor{10},rmse{10}] = mod_error(COUP050_diffDM{2},COUP050_diffDO{2});
STD{10} = nanstd(COUP050_diffDM{2});
str = {sprintf('RMSE=%.4f',round(rmse{10}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{10}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{10}.CCall(101),4))};
dim = [0.339485179407176,0.416488993339163,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,11)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (1.2>f/f_p\geq2.0)')
xlabel('Obs diffDir (^o)')
ylabel('COUP050 diffDir (^o)')
[~,Sk{11},maxcor{11},rmse{11}] = mod_error(COUP050_diffDM{3},COUP050_diffDO{3});
STD{11} = nanstd(COUP050_diffDM{3});
str = {sprintf('RMSE=%.4f',round(rmse{11}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{11}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{11}.CCall(101),4))};
dim = [0.547753510140407,0.416488993339163,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,12)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (2.0>f/f_p\geq3.0)')
xlabel('Obs diffDir (^o)')
ylabel('COUP050 diffDir (^o)')
[~,Sk{12},maxcor{12},rmse{12}] = mod_error(COUP050_diffDM{4},COUP050_diffDO{4});
STD{12} = nanstd(COUP050_diffDM{4}); 
str = {sprintf('RMSE=%.4f',round(rmse{12}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{12}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{12}.CCall(101),4))};
dim = [0.748875819032762,0.415373348783314,0.0698,0.0603];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN075 Comparison: 
SWAN075_diffDO{1} = []; SWAN075_diffDM{1} = [];
SWAN075_diffDO{2} = []; SWAN075_diffDM{2} = [];
SWAN075_diffDO{3} = []; SWAN075_diffDM{3} = [];
SWAN075_diffDO{4} = []; SWAN075_diffDM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OdiffDir1 = Buoy.Obs.%s.EWMdiffDir1;',Bnames{i}))
    eval(sprintf('MdiffDir1 = Buoy.SWAN075.%s.EWMdiffDir1;',Bnames{i}))
    eval(sprintf('OdiffDir2 = Buoy.Obs.%s.EWMdiffDir2;',Bnames{i}))
    eval(sprintf('MdiffDir2 = Buoy.SWAN075.%s.EWMdiffDir2;',Bnames{i}))
    eval(sprintf('OdiffDir3 = Buoy.Obs.%s.EWMdiffDir3;',Bnames{i}))
    eval(sprintf('MdiffDir3 = Buoy.SWAN075.%s.EWMdiffDir3;',Bnames{i}))
    eval(sprintf('OdiffDir4 = Buoy.Obs.%s.EWMdiffDir4;',Bnames{i}))
    eval(sprintf('MdiffDir4 = Buoy.SWAN075.%s.EWMdiffDir4;',Bnames{i}))
    
    subplot(4,4,13)
    scatter(OdiffDir1,MdiffDir1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,14)
    scatter(OdiffDir2,MdiffDir2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,15)
    scatter(OdiffDir3,MdiffDir3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,16)
    scatter(OdiffDir4,MdiffDir4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN075_diffDO{1} = cat(2,SWAN075_diffDO{1},OdiffDir1);
    SWAN075_diffDO{2} = cat(2,SWAN075_diffDO{2},OdiffDir2);
    SWAN075_diffDO{3} = cat(2,SWAN075_diffDO{3},OdiffDir3);
    SWAN075_diffDO{4} = cat(2,SWAN075_diffDO{4},OdiffDir4);
    
    SWAN075_diffDM{1} = cat(2,SWAN075_diffDM{1},MdiffDir1);
    SWAN075_diffDM{2} = cat(2,SWAN075_diffDM{2},MdiffDir2);
    SWAN075_diffDM{3} = cat(2,SWAN075_diffDM{3},MdiffDir3);
    SWAN075_diffDM{4} = cat(2,SWAN075_diffDM{4},MdiffDir4);
end

subplot(4,4,13)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.5>f/f_p\geq0.8)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN075 diffDir (^o)')
%legend(BTlegendnames,'NumColumns',2,'location','north')
[~,Sk{13},maxcor{13},rmse{13}] = mod_error(SWAN075_diffDM{1},SWAN075_diffDO{1});
STD{13} = nanstd(SWAN075_diffDM{1});
str = {sprintf('RMSE=%.4f',round(rmse{13}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{13}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{13}.CCall(101),4))};
dim = [0.139521060842434,0.199581126176124,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,14)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.8>f/f_p\geq1.2)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN075 diffDir (^o)')
[~,Sk{14},maxcor{14},rmse{14}] = mod_error(SWAN075_diffDM{2},SWAN075_diffDO{2});
STD{14} = nanstd(SWAN075_diffDM{2});
str = {sprintf('RMSE=%.4f',round(rmse{14}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{14}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{14}.CCall(101),4))};
dim = [0.342605304212169,0.196326768541944,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,15)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (1.2>f/f_p\geq2.0)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN075 diffDir (^o)')
[~,Sk{15},maxcor{15},rmse{15}] = mod_error(SWAN075_diffDM{3},SWAN075_diffDO{3});
STD{15} = nanstd(SWAN075_diffDM{3});
str = {sprintf('RMSE=%.4f',round(rmse{15}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{15}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{15}.CCall(101),4))};
dim = [0.54619344773791,0.199803014196637,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,16)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (2.0>f/f_p\geq3.0)')
xlabel('Obs diffDir (^o)')
ylabel('SWAN075 diffDir (^o)')
[~,Sk{16},maxcor{16},rmse{16}] = mod_error(SWAN075_diffDM{4},SWAN075_diffDO{4});
STD{16} = nanstd(SWAN075_diffDM{4}); 
str = {sprintf('RMSE=%.4f',round(rmse{16}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{16}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{16}.CCall(101),4))};
dim = [0.754336037441498,0.196369872537659,0.0698,0.0603];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


% Rename diff Direction Statistics Variables
diffDir_rmse = rmse;
diffDir_Sk = Sk;
diffDir_maxcor = maxcor;
diffDir_STD = STD;
diffDir_OSTD = OSTD;
clear rmse Sk maxcor STD OSTD



%% Plot Difference in Directional Spread b/t Buoys

figure(6);clf; % Energy Figure
set(gcf,'position',[900,50,1000,940])


    % SWAN050 Comparison: 
SWAN_diffDSO{1} = []; SWAN_diffDSM{1} = [];
SWAN_diffDSO{2} = []; SWAN_diffDSM{2} = [];
SWAN_diffDSO{3} = []; SWAN_diffDSM{3} = [];
SWAN_diffDSO{4} = []; SWAN_diffDSM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OdiffDS1 = Buoy.Obs.%s.EWMdiffDspread1;',Bnames{i}))
    eval(sprintf('MdiffDS1 = Buoy.SWAN050.%s.EWMdiffDspread1;',Bnames{i}))
    eval(sprintf('OdiffDS2 = Buoy.Obs.%s.EWMdiffDspread2;',Bnames{i}))
    eval(sprintf('MdiffDS2 = Buoy.SWAN050.%s.EWMdiffDspread2;',Bnames{i}))
    eval(sprintf('OdiffDS3 = Buoy.Obs.%s.EWMdiffDspread3;',Bnames{i}))
    eval(sprintf('MdiffDS3 = Buoy.SWAN050.%s.EWMdiffDspread3;',Bnames{i}))
    eval(sprintf('OdiffDS4 = Buoy.Obs.%s.EWMdiffDspread4;',Bnames{i}))
    eval(sprintf('MdiffDS4 = Buoy.SWAN050.%s.EWMdiffDspread4;',Bnames{i}))
    
    subplot(441)
    scatter(OdiffDS1,MdiffDS1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(442)
    scatter(OdiffDS2,MdiffDS2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(443)
    scatter(OdiffDS3,MdiffDS3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(444)
    scatter(OdiffDS4,MdiffDS4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN_diffDSO{1} = cat(2,SWAN_diffDSO{1},OdiffDS1);
    SWAN_diffDSO{2} = cat(2,SWAN_diffDSO{2},OdiffDS2);
    SWAN_diffDSO{3} = cat(2,SWAN_diffDSO{3},OdiffDS3);
    SWAN_diffDSO{4} = cat(2,SWAN_diffDSO{4},OdiffDS4);
    
    SWAN_diffDSM{1} = cat(2,SWAN_diffDSM{1},MdiffDS1);
    SWAN_diffDSM{2} = cat(2,SWAN_diffDSM{2},MdiffDS2);
    SWAN_diffDSM{3} = cat(2,SWAN_diffDSM{3},MdiffDS3);
    SWAN_diffDSM{4} = cat(2,SWAN_diffDSM{4},MdiffDS4);
end

subplot(441)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.5>f/f_p\geq0.8)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN050 diffDspread (^o)')
legend(BTlegendnames,'fontsize',7,'NumColumns',3,'position',[0.429603927017891,0.014174061193324,0.160686427457098,0.050984934922027])
[~,Sk{1},maxcor{1},rmse{1}] = mod_error(SWAN_diffDSM{1},SWAN_diffDSO{1});
STD{1} = nanstd(SWAN_diffDSM{1});
OSTD{1} = nanstd(SWAN_diffDSO{1});
str = {sprintf('RMSE=%.4f',round(rmse{1}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{1}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{1}.CCall(101),4))};
dim = [0.131996879875195,0.86558516801854,0.047410296411856,0.052178447276941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(442)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.8>f/f_p\geq1.2)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN050 diffDspread (^o)')
[~,Sk{2},maxcor{2},rmse{2}] = mod_error(SWAN_diffDSM{2},SWAN_diffDSO{2});
STD{2} = nanstd(SWAN_diffDSM{2});
OSTD{2} = nanstd(SWAN_diffDSO{2});
str = {sprintf('RMSE=%.4f',round(rmse{2}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{2}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{2}.CCall(101),4))};
dim = [0.33792511700468,0.864426419466976,0.047410296411856,0.052178447276941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(443)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (1.2>f/f_p\geq2.0)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN050 diffDspread (^o)')
[~,Sk{3},maxcor{3},rmse{3}] = mod_error(SWAN_diffDSM{3},SWAN_diffDSO{3});
STD{3} = nanstd(SWAN_diffDSM{3});
OSTD{3} = nanstd(SWAN_diffDSO{3});
str = {sprintf('RMSE=%.4f',round(rmse{3}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{3}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{3}.CCall(101),4))};
dim = [0.546973478939159,0.865585168018541,0.047410296411856,0.052178447276941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(444)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (2.0>f/f_p\geq3.0)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN050 diffDspread (^o)')
[~,Sk{4},maxcor{4},rmse{4}] = mod_error(SWAN_diffDSM{4},SWAN_diffDSO{4});
STD{4} = nanstd(SWAN_diffDSM{4});
OSTD{4} = nanstd(SWAN_diffDSO{4});
str = {sprintf('RMSE=%.4f',round(rmse{4}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{4}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{4}.CCall(101),4))};
dim = [0.752121684867396,0.86558516801854,0.047410296411856,0.052178447276941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')



    % SWAN050sm Comparison: 
SWANsm_diffDSO{1} = []; SWANsm_diffDSM{1} = [];
SWANsm_diffDSO{2} = []; SWANsm_diffDSM{2} = [];
SWANsm_diffDSO{3} = []; SWANsm_diffDSM{3} = [];
SWANsm_diffDSO{4} = []; SWANsm_diffDSM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OdiffDS1 = Buoy.Obs.%s.EWMdiffDspread1;',Bnames{i}))
    eval(sprintf('MdiffDS1 = Buoy.SWAN050sm.%s.EWMdiffDspread1;',Bnames{i}))
    eval(sprintf('OdiffDS2 = Buoy.Obs.%s.EWMdiffDspread2;',Bnames{i}))
    eval(sprintf('MdiffDS2 = Buoy.SWAN050sm.%s.EWMdiffDspread2;',Bnames{i}))
    eval(sprintf('OdiffDS3 = Buoy.Obs.%s.EWMdiffDspread3;',Bnames{i}))
    eval(sprintf('MdiffDS3 = Buoy.SWAN050sm.%s.EWMdiffDspread3;',Bnames{i}))
    eval(sprintf('OdiffDS4 = Buoy.Obs.%s.EWMdiffDspread4;',Bnames{i}))
    eval(sprintf('MdiffDS4 = Buoy.SWAN050sm.%s.EWMdiffDspread4;',Bnames{i}))
    
    subplot(445)
    scatter(OdiffDS1,MdiffDS1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(446)
    scatter(OdiffDS2,MdiffDS2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(447)
    scatter(OdiffDS3,MdiffDS3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(448)
    scatter(OdiffDS4,MdiffDS4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWANsm_diffDSO{1} = cat(2,SWANsm_diffDSO{1},OdiffDS1);
    SWANsm_diffDSO{2} = cat(2,SWANsm_diffDSO{2},OdiffDS2);
    SWANsm_diffDSO{3} = cat(2,SWANsm_diffDSO{3},OdiffDS3);
    SWANsm_diffDSO{4} = cat(2,SWANsm_diffDSO{4},OdiffDS4);
    
    SWANsm_diffDSM{1} = cat(2,SWANsm_diffDSM{1},MdiffDS1);
    SWANsm_diffDSM{2} = cat(2,SWANsm_diffDSM{2},MdiffDS2);
    SWANsm_diffDSM{3} = cat(2,SWANsm_diffDSM{3},MdiffDS3);
    SWANsm_diffDSM{4} = cat(2,SWANsm_diffDSM{4},MdiffDS4);
end

subplot(445)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.5>f/f_p\geq0.8)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN050sm diffDspread (^o)')
%legend(BTlegendnames,'NumColumns',3,'fontsize',7,'position',[0.128511883928935,0.573849611598886,0.159126363869018,0.050984934922027])
[~,Sk{5},maxcor{5},rmse{5}] = mod_error(SWANsm_diffDSM{1},SWANsm_diffDSO{1});
STD{5} = nanstd(SWANsm_diffDSM{1});
str = {sprintf('RMSE=%.4f',round(rmse{5}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{5}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{5}.CCall(101),4))};
dim = [0.133556942277692,0.646581691772885,0.047410296411856,0.052178447276941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(446)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.8>f/f_p\geq1.2)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN050sm diffDspread (^o)')
[~,Sk{6},maxcor{6},rmse{6}] = mod_error(SWANsm_diffDSM{2},SWANsm_diffDSO{2});
STD{6} = nanstd(SWANsm_diffDSM{2});
str = {sprintf('RMSE=%.4f',round(rmse{6}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{6}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{6}.CCall(101),4))};
dim = [0.341201248049923,0.639316339804979,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(447)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (1.2>f/f_p\geq2.0)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN050sm diffDspread (^o)')
[~,Sk{7},maxcor{7},rmse{7}] = mod_error(SWANsm_diffDSM{3},SWANsm_diffDSO{3});
STD{7} = nanstd(SWANsm_diffDSM{3});
str = {sprintf('RMSE=%.4f',round(rmse{7}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{7}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{7}.CCall(101),4))};
dim = [0.545569422776913,0.638157591253415,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(448)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (2.0>f/f_p\geq3.0)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN050sm diffDspread (^o)')
[~,Sk{8},maxcor{8},rmse{8}] = mod_error(SWANsm_diffDSM{4},SWANsm_diffDSO{4});
STD{8} = nanstd(SWANsm_diffDSM{4});
str = {sprintf('RMSE=%.4f',round(rmse{8}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{8}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{8}.CCall(101),4))};
dim = [0.750717628705149,0.641633836908108,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % COUP050 Comparison: 
COUP_diffDSO{1} = []; COUP_diffDSM{1} = [];
COUP_diffDSO{2} = []; COUP_diffDSM{2} = [];
COUP_diffDSO{3} = []; COUP_diffDSM{3} = [];
COUP_diffDSO{4} = []; COUP_diffDSM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OdiffDS1 = Buoy.Obs.%s.EWMdiffDspread1;',Bnames{i}))
    eval(sprintf('MdiffDS1 = Buoy.COUP050.%s.EWMdiffDspread1;',Bnames{i}))
    eval(sprintf('OdiffDS2 = Buoy.Obs.%s.EWMdiffDspread2;',Bnames{i}))
    eval(sprintf('MdiffDS2 = Buoy.COUP050.%s.EWMdiffDspread2;',Bnames{i}))
    eval(sprintf('OdiffDS3 = Buoy.Obs.%s.EWMdiffDspread3;',Bnames{i}))
    eval(sprintf('MdiffDS3 = Buoy.COUP050.%s.EWMdiffDspread3;',Bnames{i}))
    eval(sprintf('OdiffDS4 = Buoy.Obs.%s.EWMdiffDspread4;',Bnames{i}))
    eval(sprintf('MdiffDS4 = Buoy.COUP050.%s.EWMdiffDspread4;',Bnames{i}))
    
    subplot(449)
    scatter(OdiffDS1,MdiffDS1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,10)
    scatter(OdiffDS2,MdiffDS2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,11)
    scatter(OdiffDS3,MdiffDS3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,12)
    scatter(OdiffDS4,MdiffDS4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    COUP_diffDSO{1} = cat(2,COUP_diffDSO{1},OdiffDS1);
    COUP_diffDSO{2} = cat(2,COUP_diffDSO{2},OdiffDS2);
    COUP_diffDSO{3} = cat(2,COUP_diffDSO{3},OdiffDS3);
    COUP_diffDSO{4} = cat(2,COUP_diffDSO{4},OdiffDS4);
    
    COUP_diffDSM{1} = cat(2,COUP_diffDSM{1},MdiffDS1);
    COUP_diffDSM{2} = cat(2,COUP_diffDSM{2},MdiffDS2);
    COUP_diffDSM{3} = cat(2,COUP_diffDSM{3},MdiffDS3);
    COUP_diffDSM{4} = cat(2,COUP_diffDSM{4},MdiffDS4);
end

subplot(449)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.5>f/f_p\geq0.8)')
xlabel('Obs diffDspread (^o)')
ylabel('COUP050 diffDspread (^o)')
%legend(BTlegendnames,'fontsize',7,'NumColumns',3,'position',[0.129291915130183,0.27338611217826,0.159126363869018,0.050984934922027])
[~,Sk{9},maxcor{9},rmse{9}] = mod_error(COUP_diffDSM{1},COUP_diffDSO{1});
STD{9} = nanstd(COUP_diffDSM{1});
str = {sprintf('RMSE=%.4f',round(rmse{9}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{9}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{9}.CCall(101),4))};
dim = [0.132932917316694,0.419154115007761,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,10)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.8>f/f_p\geq1.2)')
xlabel('Obs diffDspread (^o)')
ylabel('COUP050 diffDspread (^o)')
[~,Sk{10},maxcor{10},rmse{10}] = mod_error(COUP_diffDSM{2},COUP_diffDSO{2});
STD{10} = nanstd(COUP_diffDSM{2});
str = {sprintf('RMSE=%.4f',round(rmse{10}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{10}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{10}.CCall(101),4))};
dim = [0.338861154446179,0.423789109214017,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,11)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (1.2>f/f_p\geq2.0)')
xlabel('Obs diffDspread (^o)')
ylabel('COUP050 diffDspread (^o)')
[~,Sk{11},maxcor{11},rmse{11}] = mod_error(COUP_diffDSM{3},COUP_diffDSO{3});
STD{11} = nanstd(COUP_diffDSM{3});
str = {sprintf('RMSE=%.4f',round(rmse{11}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{11}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{11}.CCall(101),4))};
dim = [0.544009360374417,0.421471612110889,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,12)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (2.0>f/f_p\geq3.0)')
xlabel('Obs diffDspread (^o)')
ylabel('COUP050 diffDspread (^o)')
[~,Sk{12},maxcor{12},rmse{12}] = mod_error(COUP_diffDSM{4},COUP_diffDSO{4});
STD{12} = nanstd(COUP_diffDSM{4});
str = {sprintf('RMSE=%.4f',round(rmse{12}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{12}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{12}.CCall(101),4))};
dim = [0.753837753510144,0.417995366456196,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN075 Comparison: 
SWAN075_diffDSO{1} = []; SWAN075_diffDSM{1} = [];
SWAN075_diffDSO{2} = []; SWAN075_diffDSM{2} = [];
SWAN075_diffDSO{3} = []; SWAN075_diffDSM{3} = [];
SWAN075_diffDSO{4} = []; SWAN075_diffDSM{4} = [];

for i = [1:4 6:8]
    eval(sprintf('OdiffDS1 = Buoy.Obs.%s.EWMdiffDspread1;',Bnames{i}))
    eval(sprintf('MdiffDS1 = Buoy.SWAN075.%s.EWMdiffDspread1;',Bnames{i}))
    eval(sprintf('OdiffDS2 = Buoy.Obs.%s.EWMdiffDspread2;',Bnames{i}))
    eval(sprintf('MdiffDS2 = Buoy.SWAN075.%s.EWMdiffDspread2;',Bnames{i}))
    eval(sprintf('OdiffDS3 = Buoy.Obs.%s.EWMdiffDspread3;',Bnames{i}))
    eval(sprintf('MdiffDS3 = Buoy.SWAN075.%s.EWMdiffDspread3;',Bnames{i}))
    eval(sprintf('OdiffDS4 = Buoy.Obs.%s.EWMdiffDspread4;',Bnames{i}))
    eval(sprintf('MdiffDS4 = Buoy.SWAN075.%s.EWMdiffDspread4;',Bnames{i}))
    
    subplot(4,4,13)
    scatter(OdiffDS1,MdiffDS1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,14)
    scatter(OdiffDS2,MdiffDS2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,15)
    scatter(OdiffDS3,MdiffDS3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,16)
    scatter(OdiffDS4,MdiffDS4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN075_diffDSO{1} = cat(2,SWAN075_diffDSO{1},OdiffDS1);
    SWAN075_diffDSO{2} = cat(2,SWAN075_diffDSO{2},OdiffDS2);
    SWAN075_diffDSO{3} = cat(2,SWAN075_diffDSO{3},OdiffDS3);
    SWAN075_diffDSO{4} = cat(2,SWAN075_diffDSO{4},OdiffDS4);
    
    SWAN075_diffDSM{1} = cat(2,SWAN075_diffDSM{1},MdiffDS1);
    SWAN075_diffDSM{2} = cat(2,SWAN075_diffDSM{2},MdiffDS2);
    SWAN075_diffDSM{3} = cat(2,SWAN075_diffDSM{3},MdiffDS3);
    SWAN075_diffDSM{4} = cat(2,SWAN075_diffDSM{4},MdiffDS4);
end

subplot(4,4,13)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.5>f/f_p\geq0.8)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN075 diffDspread (^o)')
%legend(BTlegendnames,'fontsize',7,'NumColumns',3,'position',[0.129291915130183,0.27338611217826,0.159126363869018,0.050984934922027])
[~,Sk{13},maxcor{13},rmse{13}] = mod_error(SWAN075_diffDSM{1},SWAN075_diffDSO{1});
STD{13} = nanstd(SWAN075_diffDSM{1});
str = {sprintf('RMSE=%.4f',round(rmse{13}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{13}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{13}.CCall(101),4))};
dim = [0.133712948517942,0.202468135865235,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,14)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.8>f/f_p\geq1.2)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN075 diffDspread (^o)')
[~,Sk{14},maxcor{14},rmse{14}] = mod_error(SWAN075_diffDSM{2},SWAN075_diffDSO{2});
STD{14} = nanstd(SWAN075_diffDSM{2});
str = {sprintf('RMSE=%.4f',round(rmse{14}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{14}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{14}.CCall(101),4))};
dim = [0.340421216848675,0.20130938731367,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,15)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (1.2>f/f_p\geq2.0)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN075 diffDspread (^o)')
[~,Sk{15},maxcor{15},rmse{15}] = mod_error(SWAN075_diffDSM{3},SWAN075_diffDSO{3});
STD{15} = nanstd(SWAN075_diffDSM{3});
str = {sprintf('RMSE=%.4f',round(rmse{15}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{15}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{15}.CCall(101),4))};
dim = [0.545569422776914,0.202468135865234,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,16)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (2.0>f/f_p\geq3.0)')
xlabel('Obs diffDspread (^o)')
ylabel('SWAN075 diffDspread (^o)')
[~,Sk{16},maxcor{16},rmse{16}] = mod_error(SWAN075_diffDSM{4},SWAN075_diffDSO{4});
STD{16} = nanstd(SWAN075_diffDSM{4});
str = {sprintf('RMSE=%.4f',round(rmse{16}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{16}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{16}.CCall(101),4))};
dim = [0.750717628705151,0.198991890210541,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')



% Rename diff Directional Spread Statistics Variables
diffDspread_rmse = rmse;
diffDspread_Sk = Sk;
diffDspread_maxcor = maxcor;
diffDspread_STD = STD;
diffDspread_OSTD = OSTD;
clear rmse Sk maxcor STD OSTD


%% Calculate Statistics for Obs and Models for each of the buoys

% Preallocate cells (to make all same size for addtion):
TEDi_Sk1 = zeros(16,7);TEDi_Sk2 = zeros(16,7);TEDi_Sk3 = zeros(16,7);TEDi_Sk4 = zeros(16,7);
TEDi_maxcor1 = zeros(16,7);TEDi_maxcor2 = zeros(16,7);TEDi_maxcor3 = zeros(16,7);TEDi_maxcor4 = zeros(16,7);
TEDi_rmse1 = zeros(16,7);TEDi_rmse2 = zeros(16,7);TEDi_rmse3 = zeros(16,7);TEDi_rmse4 = zeros(16,7);
TEDistd_M1 = zeros(16,7);TEDistd_M2 = zeros(16,7);TEDistd_M3 = zeros(16,7);TEDistd_M4 = zeros(16,7);
TEDistd_O1 = zeros(16,7);TEDistd_O2 = zeros(16,7);TEDistd_O3 = zeros(16,7);TEDistd_O4 = zeros(16,7);
dDi_Sk1 = zeros(16,7);dDi_Sk2 = zeros(16,7);dDi_Sk3 = zeros(16,7);dDi_Sk4 = zeros(16,7);
dDi_maxcor1 = zeros(16,7);dDi_maxcor2 = zeros(16,7);dDi_maxcor3 = zeros(16,7);dDi_maxcor4 = zeros(16,7);
dDi_rmse1 = zeros(16,7);dDi_rmse2 = zeros(16,7);dDi_rmse3 = zeros(16,7);dDi_rmse4 = zeros(16,7);
dDistd_M1 = zeros(16,7);dDistd_M2 = zeros(16,7);dDistd_M3 = zeros(16,7);dDistd_M4 = zeros(16,7);
dDistd_O1 = zeros(16,7);dDistd_O2 = zeros(16,7);dDistd_O3 = zeros(16,7);dDistd_O4 = zeros(16,7);
dDSi_Sk1 = zeros(16,7);dDSi_Sk2 = zeros(16,7);dDSi_Sk3 = zeros(16,7);dDSi_Sk4 = zeros(16,7);
dDSi_maxcor1 = zeros(16,7);dDSi_maxcor2 = zeros(16,7);dDSi_maxcor3 = zeros(16,7);dDSi_maxcor4 = zeros(16,7);
dDSi_rmse1 = zeros(16,7);dDSi_rmse2 = zeros(16,7);dDSi_rmse3 = zeros(16,7);dDSi_rmse4 = zeros(16,7);
dDSistd_M1 = zeros(16,7);dDSistd_M2 = zeros(16,7);dDSistd_M3 = zeros(16,7);dDSistd_M4 = zeros(16,7);
dDSistd_O1 = zeros(16,7);dDSistd_O2 = zeros(16,7);dDSistd_O3 = zeros(16,7);dDSistd_O4 = zeros(16,7);



for MO = 1:4 % for the 4 models
    
    for i = [1:4 6:8] % for b/t 9 buoys (8 total)
        
        if MO == 1
            eval(sprintf('Mod = Buoy.SWAN050.%s;',Bnames{i})) % SWAN050
        elseif MO == 2
            eval(sprintf('Mod = Buoy.SWAN050sm.%s;',Bnames{i})) % SWAN050sm
        elseif MO == 3
            eval(sprintf('Mod = Buoy.COUP050.%s;',Bnames{i})) % COUP050
        else
            eval(sprintf('Mod = Buoy.SWAN075.%s;',Bnames{i})) % SWAN075
        end
        
        eval(sprintf('Obs = Buoy.Obs.%s;',Bnames{i})) % observed
        
        % so we don't skip a line in the matrix created below
        if i>4
            j = i-1;
        else
            j = i;
        end

        % TED (4 freq bands)
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWM_TED1,Obs.EWM_TED1);
        TEDi_Sk1(MO*4-3,j) = Sk.sk1;
        TEDi_maxcor1(MO*4-3,j) = maxcor.CCall(101);
        TEDi_rmse1(MO*4-3,j) = rmse.r1;
        TEDistd_M1(MO*4-3,j) = nanstd(Mod.EWM_TED1);
        TEDistd_O1(MO*4-3,j) = nanstd(Obs.EWM_TED1);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWM_TED2,Obs.EWM_TED2);
        TEDi_Sk2(MO*4-2,j) = Sk.sk1;
        TEDi_maxcor2(MO*4-2,j) = maxcor.CCall(101);
        TEDi_rmse2(MO*4-2,j) = rmse.r1;
        TEDistd_M2(MO*4-2,j) = nanstd(Mod.EWM_TED2);
        TEDistd_O2(MO*4-2,j) = nanstd(Obs.EWM_TED2);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWM_TED3,Obs.EWM_TED3);
        TEDi_Sk3(MO*4-1,j) = Sk.sk1;
        TEDi_maxcor3(MO*4-1,j) = maxcor.CCall(101);
        TEDi_rmse3(MO*4-1,j) = rmse.r1;
        TEDistd_M3(MO*4-1,j) = nanstd(Mod.EWM_TED3);
        TEDistd_O3(MO*4-1,j) = nanstd(Obs.EWM_TED3);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWM_TED4,Obs.EWM_TED4);
        TEDi_Sk4(MO*4,j) = Sk.sk1;
        TEDi_maxcor4(MO*4,j) = maxcor.CCall(101);
        TEDi_rmse4(MO*4,j) = rmse.r1;
        TEDistd_M4(MO*4,j) = nanstd(Mod.EWM_TED4);
        TEDistd_O4(MO*4,j) = nanstd(Obs.EWM_TED4);
        
        
        % Diff Direction (4 freq bands)
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMdiffDir1,Obs.EWMdiffDir1);
        dDi_Sk1(MO*4-3,j) = Sk.sk1;
        dDi_maxcor1(MO*4-3,j) = maxcor.CCall(101);
        dDi_rmse1(MO*4-3,j) = rmse.r1;
        dDistd_M1(MO*4-3,j) = nanstd(Mod.EWMdiffDir1);
        dDistd_O1(MO*4-3,j) = nanstd(Obs.EWMdiffDir1);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMdiffDir2,Obs.EWMdiffDir2);
        dDi_Sk2(MO*4-2,j) = Sk.sk1;
        dDi_maxcor2(MO*4-2,j) = maxcor.CCall(101);
        dDi_rmse2(MO*4-2,j) = rmse.r1;
        dDistd_M2(MO*4-2,j) = nanstd(Mod.EWMdiffDir2);
        dDistd_O2(MO*4-2,j) = nanstd(Obs.EWMdiffDir2);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMdiffDir3,Obs.EWMdiffDir3);
        dDi_Sk3(MO*4-1,j) = Sk.sk1;
        dDi_maxcor3(MO*4-1,j) = maxcor.CCall(101);
        dDi_rmse3(MO*4-1,j) = rmse.r1;
        dDistd_M3(MO*4-1,j) = nanstd(Mod.EWMdiffDir3);
        dDistd_O3(MO*4-1,j) = nanstd(Obs.EWMdiffDir3);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMdiffDir4,Obs.EWMdiffDir4);
        dDi_Sk4(MO*4,j) = Sk.sk1;
        dDi_maxcor4(MO*4,j) = maxcor.CCall(101);
        dDi_rmse4(MO*4,j) = rmse.r1;
        dDistd_M4(MO*4,j) = nanstd(Mod.EWMdiffDir4);
        dDistd_O4(MO*4,j) = nanstd(Obs.EWMdiffDir4);
        
        
        % Diff Directional Spread (4 freq bands)
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMdiffDspread1,Obs.EWMdiffDspread1);
        dDSi_Sk1(MO*4-3,j) = Sk.sk1;
        dDSi_maxcor1(MO*4-3,j) = maxcor.CCall(101);
        dDSi_rmse1(MO*4-3,j) = rmse.r1;
        dDSistd_M1(MO*4-3,j) = nanstd(Mod.EWMdiffDspread1);
        dDSistd_O1(MO*4-3,j) = nanstd(Obs.EWMdiffDspread1);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMdiffDspread2,Obs.EWMdiffDspread2);
        dDSi_Sk2(MO*4-2,j) = Sk.sk1;
        dDSi_maxcor2(MO*4-2,j) = maxcor.CCall(101);
        dDSi_rmse2(MO*4-2,j) = rmse.r1;
        dDSistd_M2(MO*4-2,j) = nanstd(Mod.EWMdiffDspread2);
        dDSistd_O2(MO*4-2,j) = nanstd(Obs.EWMdiffDspread2);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMdiffDspread3,Obs.EWMdiffDspread3);
        dDSi_Sk3(MO*4-1,j) = Sk.sk1;
        dDSi_maxcor3(MO*4-1,j) = maxcor.CCall(101);
        dDSi_rmse3(MO*4-1,j) = rmse.r1;
        dDSistd_M3(MO*4-1,j) = nanstd(Mod.EWMdiffDspread3);
        dDSistd_O3(MO*4-1,j) = nanstd(Obs.EWMdiffDspread3);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMdiffDspread4,Obs.EWMdiffDspread4);
        dDSi_Sk4(MO*4,j) = Sk.sk1;
        dDSi_maxcor4(MO*4,j) = maxcor.CCall(101);
        dDSi_rmse4(MO*4,j) = rmse.r1;
        dDSistd_M4(MO*4,j) = nanstd(Mod.EWMdiffDspread4);
        dDSistd_O4(MO*4,j) = nanstd(Obs.EWMdiffDspread4);

    end
end

% Add all together

TEDi_Sk = TEDi_Sk1 + TEDi_Sk2 + TEDi_Sk3 + TEDi_Sk4;
TEDi_maxcor = TEDi_maxcor1 + TEDi_maxcor2 + TEDi_maxcor3 + TEDi_maxcor4;
TEDi_rmse = TEDi_rmse1 + TEDi_rmse2 + TEDi_rmse3 + TEDi_rmse4;
TEDi_std_M = TEDistd_M1 + TEDistd_M2 + TEDistd_M3 + TEDistd_M4;
TEDi_std_O = TEDistd_O1 + TEDistd_O2 + TEDistd_O3 + TEDistd_O4;
dDi_Sk = dDi_Sk1 + dDi_Sk2 + dDi_Sk3 + dDi_Sk4;
dDi_maxcor = dDi_maxcor1 + dDi_maxcor2 + dDi_maxcor3 + dDi_maxcor4;
dDi_rmse = dDi_rmse1 + dDi_rmse2 + dDi_rmse3 + dDi_rmse4;
dDi_std_M = dDistd_M1 + dDistd_M2 + dDistd_M3 + dDistd_M4;
dDi_std_O = dDistd_O1 + dDistd_O2 + dDistd_O3 + dDistd_O4;
dDSi_Sk = dDSi_Sk1 + dDSi_Sk2 + dDSi_Sk3 + dDSi_Sk4;
dDSi_maxcor = dDSi_maxcor1 + dDSi_maxcor2 + dDSi_maxcor3 + dDSi_maxcor4;
dDSi_rmse = dDSi_rmse1 + dDSi_rmse2 + dDSi_rmse3 + dDSi_rmse4;
dDSi_std_M = dDSistd_M1 + dDSistd_M2 + dDSistd_M3 + dDSistd_M4;
dDSi_std_O = dDSistd_O1 + dDSistd_O2 + dDSistd_O3 + dDSistd_O4;

% Clear up individual frequency variables
clear TEDi_Sk1 + TEDi_Sk2 + TEDi_Sk3 + TEDi_Sk4 TEDi_maxcor1 + TEDi_maxcor2 + TEDi_maxcor3 + TEDi_maxcor4
clear TEDi_rmse1 + TEDi_rmse2 + TEDi_rmse3 + TEDi_rmse4 TEDistd_M1 + TEDistd_M2 + TEDistd_M3 + TEDistd_M4
clear TEDistd_O1 + TEDistd_O2 + TEDistd_O3 + TEDistd_O4 dDi_Sk1 + dDi_Sk2 + dDi_Sk3 + dDi_Sk4
clear dDi_maxcor1 + dDi_maxcor2 + dDi_maxcor3 + dDi_maxcor4 dDi_rmse1 + dDi_rmse2 + dDi_rmse3 + dDi_rmse4
clear dDistd_M1 + dDistd_M2 + dDistd_M3 + dDistd_M4 dDistd_O1 + dDistd_O2 + dDistd_O3 + dDistd_O4
clear dDSi_Sk1 + dDSi_Sk2 + dDSi_Sk3 + dDSi_Sk4 dDSi_maxcor1 + dDSi_maxcor2 + dDSi_maxcor3 + dDSi_maxcor4
clear dDSi_rmse1 + dDSi_rmse2 + dDSi_rmse3 + dDSi_rmse4 dDSistd_M1 + dDSistd_M2 + dDSistd_M3 + dDSistd_M4
clear dDSistd_O1 + dDSistd_O2 + dDSistd_O3 + dDSistd_O4



%% Statistics Tables

        
Model = ["SWAN050";"SWAN050";"SWAN050";"SWAN050";"SWAN050sm";"SWAN050sm";...
    "SWAN050sm";"SWAN050sm";"COUP050";"COUP050";"COUP050";"COUP050";...
    "SWAN075";"SWAN075";"SWAN075";"SWAN075"];
f_Band = ["f/fp=0.5-0.8";"f/fp=0.8-1.2";"f/fp=1.2-2.0";"f/fp=2.0-3.0";...
    "f/fp=0.5-0.8";"f/fp=0.8-1.2";"f/fp=1.2-2.0";"f/fp=2.0-3.0";...
    "f/fp=0.5-0.8";"f/fp=0.8-1.2";"f/fp=1.2-2.0";"f/fp=2.0-3.0";...
    "f/fp=0.5-0.8";"f/fp=0.8-1.2";"f/fp=1.2-2.0";"f/fp=2.0-3.0"];



    % TED TABLE:
RMSE = [TED_rmse{1}.r1;TED_rmse{2}.r1;TED_rmse{3}.r1;TED_rmse{4}.r1;TED_rmse{5}.r1;...
    TED_rmse{6}.r1;TED_rmse{7}.r1;TED_rmse{8}.r1;TED_rmse{9}.r1;TED_rmse{10}.r1;...
    TED_rmse{11}.r1;TED_rmse{12}.r1;TED_rmse{13}.r1;TED_rmse{14}.r1;...
    TED_rmse{15}.r1;TED_rmse{16}.r1];
RMSE = cat(2,RMSE,TEDi_rmse);
RMSE = round(RMSE,3);

SS = [TED_Sk{1}.sk1;TED_Sk{2}.sk1;TED_Sk{3}.sk1;TED_Sk{4}.sk1;TED_Sk{5}.sk1;TED_Sk{6}.sk1;...
    TED_Sk{7}.sk1;TED_Sk{8}.sk1;TED_Sk{9}.sk1;TED_Sk{10}.sk1;TED_Sk{11}.sk1;...
    TED_Sk{12}.sk1;TED_Sk{13}.sk1;TED_Sk{14}.sk1;TED_Sk{15}.sk1;TED_Sk{16}.sk1];
SS = cat(2,SS,TEDi_Sk);
SS = round(SS,3);

MaxCor = [TED_maxcor{1}.CCall(101);TED_maxcor{2}.CCall(101);TED_maxcor{3}.CCall(101);...
    TED_maxcor{4}.CCall(101);TED_maxcor{5}.CCall(101);TED_maxcor{6}.CCall(101);...
    TED_maxcor{7}.CCall(101);TED_maxcor{8}.CCall(101);TED_maxcor{9}.CCall(101);...
    TED_maxcor{10}.CCall(101);TED_maxcor{11}.CCall(101);TED_maxcor{12}.CCall(101);
    TED_maxcor{13}.CCall(101);TED_maxcor{14}.CCall(101);TED_maxcor{15}.CCall(101);...
    TED_maxcor{16}.CCall(101)];
MaxCor = cat(2,MaxCor,TEDi_maxcor);
MaxCor = round(MaxCor,3);

STD = [TED_STD{1};TED_STD{2};TED_STD{3};TED_STD{4};TED_STD{5};TED_STD{6};...
    TED_STD{7};TED_STD{8};TED_STD{9};TED_STD{10};TED_STD{11};TED_STD{12};...
    TED_STD{13};TED_STD{14};TED_STD{15};TED_STD{16}];
STD = cat(2,STD,TEDi_std_M);
STD = round(STD,3);

TEDStats = table(Model,f_Band,RMSE,SS,MaxCor,STD)

clear RMSE SS MaxCor STD


    % diffDir TABLE:
RMSE = [diffDir_rmse{1}.r1;diffDir_rmse{2}.r1;diffDir_rmse{3}.r1;diffDir_rmse{4}.r1;diffDir_rmse{5}.r1;...
    diffDir_rmse{6}.r1;diffDir_rmse{7}.r1;diffDir_rmse{8}.r1;diffDir_rmse{9}.r1;diffDir_rmse{10}.r1;...
    diffDir_rmse{11}.r1;diffDir_rmse{12}.r1;diffDir_rmse{13}.r1;...
    diffDir_rmse{14}.r1;diffDir_rmse{15}.r1;diffDir_rmse{16}.r1];
RMSE = cat(2,RMSE,dDi_rmse);
RMSE = round(RMSE,3);

SS = [diffDir_Sk{1}.sk1;diffDir_Sk{2}.sk1;diffDir_Sk{3}.sk1;diffDir_Sk{4}.sk1;diffDir_Sk{5}.sk1;diffDir_Sk{6}.sk1;...
    diffDir_Sk{7}.sk1;diffDir_Sk{8}.sk1;diffDir_Sk{9}.sk1;diffDir_Sk{10}.sk1;diffDir_Sk{11}.sk1;...
    diffDir_Sk{12}.sk1;diffDir_Sk{13}.sk1;diffDir_Sk{14}.sk1;diffDir_Sk{15}.sk1;diffDir_Sk{16}.sk1];
SS = cat(2,SS,dDi_Sk);
SS = round(SS,3);

MaxCor = [diffDir_maxcor{1}.CCall(101);diffDir_maxcor{2}.CCall(101);diffDir_maxcor{3}.CCall(101);...
    diffDir_maxcor{4}.CCall(101);diffDir_maxcor{5}.CCall(101);diffDir_maxcor{6}.CCall(101);...
    diffDir_maxcor{7}.CCall(101);diffDir_maxcor{8}.CCall(101);diffDir_maxcor{9}.CCall(101);...
    diffDir_maxcor{10}.CCall(101);diffDir_maxcor{11}.CCall(101);diffDir_maxcor{12}.CCall(101);...
    diffDir_maxcor{13}.CCall(101);diffDir_maxcor{14}.CCall(101);diffDir_maxcor{15}.CCall(101);...
    diffDir_maxcor{16}.CCall(101)];
MaxCor = cat(2,MaxCor,dDi_maxcor);
MaxCor = round(MaxCor,3);

STD = [diffDir_STD{1};diffDir_STD{2};diffDir_STD{3};diffDir_STD{4};diffDir_STD{5};...
    diffDir_STD{6};diffDir_STD{7};diffDir_STD{8};diffDir_STD{9};...
    diffDir_STD{10};diffDir_STD{11};diffDir_STD{12};diffDir_STD{13};...
    diffDir_STD{14};diffDir_STD{15};diffDir_STD{16}];
STD = cat(2,STD,dDi_std_M);
STD = round(STD,3);

diffDirStats = table(Model,f_Band,RMSE,SS,MaxCor,STD)

clear RMSE SS MaxCor STD


    % diffDspread TABLE:
RMSE = [diffDspread_rmse{1}.r1;diffDspread_rmse{2}.r1;diffDspread_rmse{3}.r1;diffDspread_rmse{4}.r1;diffDspread_rmse{5}.r1;...
    diffDspread_rmse{6}.r1;diffDspread_rmse{7}.r1;diffDspread_rmse{8}.r1;diffDspread_rmse{9}.r1;diffDspread_rmse{10}.r1;...
    diffDspread_rmse{11}.r1;diffDspread_rmse{12}.r1;diffDspread_rmse{13}.r1;diffDspread_rmse{14}.r1;...
    diffDspread_rmse{15}.r1;diffDspread_rmse{16}.r1];
RMSE = cat(2,RMSE,dDSi_rmse);
RMSE = round(RMSE,3);

SS = [diffDspread_Sk{1}.sk1;diffDspread_Sk{2}.sk1;diffDspread_Sk{3}.sk1;diffDspread_Sk{4}.sk1;diffDspread_Sk{5}.sk1;diffDspread_Sk{6}.sk1;...
    diffDspread_Sk{7}.sk1;diffDspread_Sk{8}.sk1;diffDspread_Sk{9}.sk1;diffDspread_Sk{10}.sk1;diffDspread_Sk{11}.sk1;...
    diffDspread_Sk{12}.sk1;diffDspread_Sk{13}.sk1;diffDspread_Sk{14}.sk1;diffDspread_Sk{15}.sk1;diffDspread_Sk{16}.sk1];
SS = cat(2,SS,dDSi_Sk);
SS = round(SS,3);

MaxCor = [diffDspread_maxcor{1}.CCall(101);diffDspread_maxcor{2}.CCall(101);diffDspread_maxcor{3}.CCall(101);...
    diffDspread_maxcor{4}.CCall(101);diffDspread_maxcor{5}.CCall(101);diffDspread_maxcor{6}.CCall(101);...
    diffDspread_maxcor{7}.CCall(101);diffDspread_maxcor{8}.CCall(101);diffDspread_maxcor{9}.CCall(101);...
    diffDspread_maxcor{10}.CCall(101);diffDspread_maxcor{11}.CCall(101);diffDspread_maxcor{12}.CCall(101);...
    diffDspread_maxcor{13}.CCall(101);diffDspread_maxcor{14}.CCall(101);diffDspread_maxcor{15}.CCall(101);...
    diffDspread_maxcor{16}.CCall(101)];
MaxCor = cat(2,MaxCor,dDSi_maxcor);
MaxCor = round(MaxCor,3);

STD = [diffDspread_STD{1};diffDspread_STD{2};diffDspread_STD{3};diffDspread_STD{4};diffDspread_STD{5};...
    diffDspread_STD{6};diffDspread_STD{7};diffDspread_STD{8};diffDspread_STD{9};...
    diffDspread_STD{10};diffDspread_STD{11};diffDspread_STD{12};diffDspread_STD{13};...
    diffDspread_STD{14};diffDspread_STD{15};diffDspread_STD{16}];
STD = cat(2,STD,dDSi_std_M);
STD = round(STD,3);

diffDspreadStats = table(Model,f_Band,RMSE,SS,MaxCor,STD)

clear RMSE SS MaxCor STD


% OBSERVED STD FOR INDIVIDUAL BUOYS

TEDi_STD = TEDi_std_O;
dDi_STD = dDi_std_O;
dDSi_STD = dDSi_std_O;

iObsSTD_BT = table(Model,f_Band,TEDi_STD,dDi_STD,dDSi_STD)

clear TEDi_STD dDi_STD dDSi_STD


%% Save and Clear Variables

save('InterpModelObs.mat','-append','Buoy')

save('obsVSmod_STATS.mat','-append','TEDStats','TED_OSTD',...
    'diffDirStats','diffDir_OSTD','diffDspreadStats','diffDspread_OSTD',...
    'iObsSTD_BT')

saveas(figure(4),'Figures/ObsvsMod_TED.jpeg')
saveas(figure(5),'Figures/ObsvsMod_diffDir.jpeg')
saveas(figure(6),'Figures/ObsvsMod_diffDspread.jpeg')


clearvars -except keepVariables Buoy TEDStats iObsSTD_BT diffDspreadStats diffDirStats 







