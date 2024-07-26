%% plot_MvsO.m
%
% Noah Clark        5/17/24
%
% Purpose: -To integrate the energy per hour over the 4 frequency bands for
%          the observed data and for the 3 models. Then plot (12 plots).
%          -To determine the energy weighted mean directions and spreads for
%          the models and for the observed in each of the 4 frequency 
%          bands. Then plot each of the modeled directions and spreads
%          against the observed.
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%
% TO DO :
%       - finalize plots
%       - make sure the first plot looks right for 1st and 4th columns (WHY DO THEY NOT???)
%       - add new model to the BT version of script and then add to taylor
%         diagrams


%% Preliminaries

calc_pFreq % this script needs to be ran first to calculate the peak frequencies and update InterpModelObs.mat

clc;clear;

load('InterpModelObs.mat')

addpath('../Stat Functions') % the functions for the statistical analysis
% the function for TED analysis:
addpath('C:\Users\nsc4750\OneDrive - UNC-Wilmington\CMS Work\Summer (0516-0728)\Main Functions') 


%% Calculate
% - Do all of the integrating and energy weighted means per frequency band
%   and save everything to varibales to be plotted later on.


Bnames = fields(Buoy.Obs); % to list buoy names


%----------------------- FOR E, DIR, & DSPREAD ---------------------------
for MO = 1:5 % for observed and 4 models
    
    for i = 1:9 % for the 9 buoys
        
        if MO == 1
            eval(sprintf('STODO = Buoy.Obs.%s;',Bnames{i})) % Observed
        elseif MO == 2
            eval(sprintf('STODO = Buoy.SWAN050.%s;',Bnames{i})) % SWAN050
        elseif MO == 3
            eval(sprintf('STODO = Buoy.SWAN050sm.%s;',Bnames{i})) % SWAN050sm
        elseif MO == 4
            eval(sprintf('STODO = Buoy.COUP050.%s;',Bnames{i})) % COUP050
        else 
            eval(sprintf('STODO = Buoy.SWAN075.%s;',Bnames{i})) % SWAN075
        end
        
        
        for j = 1:length(Buoy.time) % for each individual hour
            fdivfp = Buoy.freq./STODO.pfreq(j);
            
            ind = find(fdivfp>0.5 & fdivfp<=0.8);
            intE1(j) = sum(STODO.See(ind,j)).*0.01; % integrating See
            m0 = sum(STODO.See(ind,j)).*0.01; % m0
            Dm1 = sum(STODO.See(ind,j).*STODO.Dir(ind,j)).*0.01; % Dir m1
            DSm1 = sum(STODO.See(ind,j).*STODO.Dspread(ind,j)).*0.01; %Dir spread m1
            D_EWM1(j) = Dm1/m0; % EWM direction
            DS_EWM1(j) = DSm1/m0; % EWM directional spread
            
            ind = find(fdivfp>0.8 & fdivfp<=1.2);
            intE2(j) = sum(STODO.See(ind,j)).*0.01;
            m0 = sum(STODO.See(ind,j)).*0.01;
            Dm1 = sum(STODO.See(ind,j).*STODO.Dir(ind,j)).*0.01;
            DSm1 = sum(STODO.See(ind,j).*STODO.Dspread(ind,j)).*0.01; 
            D_EWM2(j) = Dm1/m0; 
            DS_EWM2(j) = DSm1/m0;
            
            ind = find(fdivfp>1.2 & fdivfp<=2.0);
            intE3(j) = sum(STODO.See(ind,j)).*0.01;
            m0 = sum(STODO.See(ind,j)).*0.01;
            Dm1 = sum(STODO.See(ind,j).*STODO.Dir(ind,j)).*0.01;
            DSm1 = sum(STODO.See(ind,j).*STODO.Dspread(ind,j)).*0.01; 
            D_EWM3(j) = Dm1/m0; 
            DS_EWM3(j) = DSm1/m0;
            
            ind = find(fdivfp>2.0 & fdivfp<=3.0);
            intE4(j) = sum(STODO.See(ind,j)).*0.01;
            m0 = sum(STODO.See(ind,j)).*0.01;
            Dm1 = sum(STODO.See(ind,j).*STODO.Dir(ind,j)).*0.01;
            DSm1 = sum(STODO.See(ind,j).*STODO.Dspread(ind,j)).*0.01; 
            D_EWM4(j) = Dm1/m0; 
            DS_EWM4(j) = DSm1/m0;
        end
        
        if MO == 1
            eval(sprintf('Buoy.Obs.%s.intE1 = intE1;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.intE2 = intE2;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.intE3 = intE3;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.intE4 = intE4;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMDir1 = D_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMDir2 = D_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMDir3 = D_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMDir4 = D_EWM4;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMDspread1 = DS_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMDspread2 = DS_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMDspread3 = DS_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMDspread4 = DS_EWM4;',Bnames{i}))
        elseif MO == 2
            eval(sprintf('Buoy.SWAN050.%s.intE1 = intE1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.intE2 = intE2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.intE3 = intE3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.intE4 = intE4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMDir1 = D_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMDir2 = D_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMDir3 = D_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMDir4 = D_EWM4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMDspread1 = DS_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMDspread2 = DS_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMDspread3 = DS_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMDspread4 = DS_EWM4;',Bnames{i}))
        elseif MO == 3
            eval(sprintf('Buoy.SWAN050sm.%s.intE1 = intE1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.intE2 = intE2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.intE3 = intE3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.intE4 = intE4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMDir1 = D_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMDir2 = D_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMDir3 = D_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMDir4 = D_EWM4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMDspread1 = DS_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMDspread2 = DS_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMDspread3 = DS_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMDspread4 = DS_EWM4;',Bnames{i}))
        elseif MO == 4
            eval(sprintf('Buoy.COUP050.%s.intE1 = intE1;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.intE2 = intE2;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.intE3 = intE3;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.intE4 = intE4;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMDir1 = D_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMDir2 = D_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMDir3 = D_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMDir4 = D_EWM4;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMDspread1 = DS_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMDspread2 = DS_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMDspread3 = DS_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMDspread4 = DS_EWM4;',Bnames{i}))
        else
            eval(sprintf('Buoy.SWAN075.%s.intE1 = intE1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.intE2 = intE2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.intE3 = intE3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.intE4 = intE4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMDir1 = D_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMDir2 = D_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMDir3 = D_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMDir4 = D_EWM4;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMDspread1 = DS_EWM1;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMDspread2 = DS_EWM2;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMDspread3 = DS_EWM3;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMDspread4 = DS_EWM4;',Bnames{i}))
        end
             
    end
end

clear i j MO ind intE1 intE2 intE3 intE4 STODO D_EWM1 D_EWM2 D_EWM3 D_EWM4 
clear DS_EWM1 DS_EWM2 DS_EWM3 DS_EWM4 m0 Dm1 MSm1 DSm1


%% Plot Energy Comparison

lightgrey = [0.85 0.85 0.85]; %the color of the background of subplots
pformat = {'xb','xm','xr','xc','xk','ob','om','or','oc'}; % The color and shape of markers
transformat = [1:-1/11:0.27]; %the transparency of the markers
MSize = 35; %MarkerSize

one_to_one = -30:360; % to plot a one to one line through plots

    
figure(1);clf; % Energy Figure
set(gcf,'position',[150,50,1000,940])


    % SWAN050 Comparison:
    
SWAN_CATO{1} = []; SWAN_CATM{1} = [];
SWAN_CATO{2} = []; SWAN_CATM{2} = [];
SWAN_CATO{3} = []; SWAN_CATM{3} = [];
SWAN_CATO{4} = []; SWAN_CATM{4} = [];

for i = 1:9
    eval(sprintf('OintE1 = Buoy.Obs.%s.intE1;',Bnames{i}))
    eval(sprintf('MintE1 = Buoy.SWAN050.%s.intE1;',Bnames{i}))
    eval(sprintf('OintE2 = Buoy.Obs.%s.intE2;',Bnames{i}))
    eval(sprintf('MintE2 = Buoy.SWAN050.%s.intE2;',Bnames{i}))
    eval(sprintf('OintE3 = Buoy.Obs.%s.intE3;',Bnames{i}))
    eval(sprintf('MintE3 = Buoy.SWAN050.%s.intE3;',Bnames{i}))
    eval(sprintf('OintE4 = Buoy.Obs.%s.intE4;',Bnames{i}))
    eval(sprintf('MintE4 = Buoy.SWAN050.%s.intE4;',Bnames{i}))
    
    subplot(441)
    scatter(OintE1,MintE1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(442)
    scatter(OintE2,MintE2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(443)
    scatter(OintE3,MintE3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(444)
    scatter(OintE4,MintE4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN_CATO{1} = cat(2,SWAN_CATO{1},OintE1);
    SWAN_CATO{2} = cat(2,SWAN_CATO{2},OintE2);
    SWAN_CATO{3} = cat(2,SWAN_CATO{3},OintE3);
    SWAN_CATO{4} = cat(2,SWAN_CATO{4},OintE4);
    
    SWAN_CATM{1} = cat(2,SWAN_CATM{1},MintE1);
    SWAN_CATM{2} = cat(2,SWAN_CATM{2},MintE2);
    SWAN_CATM{3} = cat(2,SWAN_CATM{3},MintE3);
    SWAN_CATM{4} = cat(2,SWAN_CATM{4},MintE4);
end
    
subplot(441)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([0 xL(2)]); ylim([0 yL(2)])
axis equal; grid on;
title('Obs vs SWAN050 (0.5>f/f_p\geq0.8)')
xlabel('Obs E (m^2)')
ylabel('SWAN050 E (m^2)')
legend(Bnames,'NumColumns',2)
[~,Sk{1},maxcor{1},rmse{1}] = mod_error(SWAN_CATM{1},SWAN_CATO{1});
STD{1} = std(SWAN_CATM{1});
OSTD{1} = std(SWAN_CATO{1});
str = {sprintf('RMSE=%.4f',round(rmse{1}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{1}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{1}.CCall(101),4))};
dim = [0.212,0.805223068552775,0.057,0.045038084874864];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

% a = gca; xL = a.XLim; yL = a.YLim;
% set(a,'color',lightgrey)
% xlim([0 xL(2)]); ylim([0 yL(2)])
% plot(one_to_one,one_to_one,'--k')
% axis equal; grid on;


subplot(442)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.8>f/f_p\geq1.2)')
xlabel('Obs E (m^2)')
ylabel('SWAN050 E (m^2)')
[~,Sk{2},maxcor{2},rmse{2}] = mod_error(SWAN_CATM{2},SWAN_CATO{2});
OSTD{2} = nanstd(SWAN_CATO{2});
STD{2} = std(SWAN_CATM{2});
str = {sprintf('RMSE=%.4f',round(rmse{2}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{2}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{2}.CCall(101),4))};
dim = [0.401151119141978,0.860944709017758,0.067445363700615,0.063882062344352];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(443)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (1.2>f/f_p\geq2.0)')
xlabel('Obs E (m^2)')
ylabel('SWAN050 E (m^2)')
[~,Sk{3},maxcor{3},rmse{3}] = mod_error(SWAN_CATM{3},SWAN_CATO{3});
STD{3} = std(SWAN_CATM{3});
OSTD{3} = nanstd(SWAN_CATO{3});
str = {sprintf('RMSE=%.4f',round(rmse{3}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{3}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{3}.CCall(101),4))};
dim = [0.607488299531983,0.858778105373347,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(444)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (2.0>f/f_p\geq3.0)')
xlabel('Obs E (m^2)')
ylabel('SWAN050 E (m^2)')
[~,Sk{4},maxcor{4},rmse{4}] = mod_error(SWAN_CATM{4},SWAN_CATO{4});
STD{4} = std(SWAN_CATM{4});
OSTD{4} = nanstd(SWAN_CATO{4});
str = {sprintf('RMSE=%.4f',round(rmse{4}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{4}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{4}.CCall(101),4))};
dim = [0.811232449297974,0.861095602476476,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN050sm Comparison:
    
SWANsm_CATO{1} = []; SWANsm_CATM{1} = [];
SWANsm_CATO{2} = []; SWANsm_CATM{2} = [];
SWANsm_CATO{3} = []; SWANsm_CATM{3} = [];
SWANsm_CATO{4} = []; SWANsm_CATM{4} = [];

for i = 1:9
    eval(sprintf('OintE1 = Buoy.Obs.%s.intE1;',Bnames{i}))
    eval(sprintf('MintE1 = Buoy.SWAN050sm.%s.intE1;',Bnames{i}))
    eval(sprintf('OintE2 = Buoy.Obs.%s.intE2;',Bnames{i}))
    eval(sprintf('MintE2 = Buoy.SWAN050sm.%s.intE2;',Bnames{i}))
    eval(sprintf('OintE3 = Buoy.Obs.%s.intE3;',Bnames{i}))
    eval(sprintf('MintE3 = Buoy.SWAN050sm.%s.intE3;',Bnames{i}))
    eval(sprintf('OintE4 = Buoy.Obs.%s.intE4;',Bnames{i}))
    eval(sprintf('MintE4 = Buoy.SWAN050sm.%s.intE4;',Bnames{i}))
    
    subplot(445)
    scatter(OintE1,MintE1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(446)
    scatter(OintE2,MintE2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(447)
    scatter(OintE3,MintE3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(448)
    scatter(OintE4,MintE4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWANsm_CATO{1} = cat(2,SWANsm_CATO{1},OintE1);
    SWANsm_CATO{2} = cat(2,SWANsm_CATO{2},OintE2);
    SWANsm_CATO{3} = cat(2,SWANsm_CATO{3},OintE3);
    SWANsm_CATO{4} = cat(2,SWANsm_CATO{4},OintE4);
    
    SWANsm_CATM{1} = cat(2,SWANsm_CATM{1},MintE1);
    SWANsm_CATM{2} = cat(2,SWANsm_CATM{2},MintE2);
    SWANsm_CATM{3} = cat(2,SWANsm_CATM{3},MintE3);
    SWANsm_CATM{4} = cat(2,SWANsm_CATM{4},MintE4);
end
    
subplot(445)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.5>f/f_p\geq0.8)')
xlabel('Obs E (m^2)')
ylabel('SWAN050sm E (m^2)')
legend(Bnames,'NumColumns',2)
[~,Sk{5},maxcor{5},rmse{5}] = mod_error(SWANsm_CATM{1},SWANsm_CATO{1});
STD{5} = std(SWANsm_CATM{1});
str = {sprintf('RMSE=%.4f',round(rmse{5}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{5}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{5}.CCall(101),4))};
dim = [0.204659906396256,0.569259335428466,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(446)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.8>f/f_p\geq1.2)')
xlabel('Obs E (m^2)')
ylabel('SWAN050sm E (m^2)')
[~,Sk{6},maxcor{6},rmse{6}] = mod_error(SWANsm_CATM{2},SWANsm_CATO{2});
STD{6} = std(SWANsm_CATM{2});
str = {sprintf('RMSE=%.4f',round(rmse{6}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{6}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{6}.CCall(101),4))};
dim = [0.404603744149766,0.645819957899613,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(447)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (1.2>f/f_p\geq2.0)')
xlabel('Obs E (m^2)')
ylabel('SWAN050sm E (m^2)')
[~,Sk{7},maxcor{7},rmse{7}] = mod_error(SWANsm_CATM{3},SWANsm_CATO{3});
STD{7} = std(SWANsm_CATM{3});
str = {sprintf('RMSE=%.4f',round(rmse{7}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{7}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{7}.CCall(101),4))};
dim = [0.610347893915758,0.643220023717539,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(448)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (2.0>f/f_p\geq3.0)')
xlabel('Obs E (m^2)')
ylabel('SWAN050sm E (m^2)')
[~,Sk{8},maxcor{8},rmse{8}] = mod_error(SWANsm_CATM{4},SWANsm_CATO{4});
STD{8} = std(SWANsm_CATM{4});
str = {sprintf('RMSE=%.4f',round(rmse{8}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{8}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{8}.CCall(101),4))};
dim = [0.807531981279252,0.640026215141791,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % COUP050 Comparison:
    
COUP_CATO{1} = []; COUP_CATM{1} = [];
COUP_CATO{2} = []; COUP_CATM{2} = [];
COUP_CATO{3} = []; COUP_CATM{3} = [];
COUP_CATO{4} = []; COUP_CATM{4} = [];

for i = 1:9
    eval(sprintf('OintE1 = Buoy.Obs.%s.intE1;',Bnames{i}))
    eval(sprintf('MintE1 = Buoy.COUP050.%s.intE1;',Bnames{i}))
    eval(sprintf('OintE2 = Buoy.Obs.%s.intE2;',Bnames{i}))
    eval(sprintf('MintE2 = Buoy.COUP050.%s.intE2;',Bnames{i}))
    eval(sprintf('OintE3 = Buoy.Obs.%s.intE3;',Bnames{i}))
    eval(sprintf('MintE3 = Buoy.COUP050.%s.intE3;',Bnames{i}))
    eval(sprintf('OintE4 = Buoy.Obs.%s.intE4;',Bnames{i}))
    eval(sprintf('MintE4 = Buoy.COUP050.%s.intE4;',Bnames{i}))
    
    subplot(449)
    scatter(OintE1,MintE1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,10)
    scatter(OintE2,MintE2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,11)
    scatter(OintE3,MintE3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,12)
    scatter(OintE4,MintE4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    COUP_CATO{1} = cat(2,COUP_CATO{1},OintE1);
    COUP_CATO{2} = cat(2,COUP_CATO{2},OintE2);
    COUP_CATO{3} = cat(2,COUP_CATO{3},OintE3);
    COUP_CATO{4} = cat(2,COUP_CATO{4},OintE4);
    
    COUP_CATM{1} = cat(2,COUP_CATM{1},MintE1);
    COUP_CATM{2} = cat(2,COUP_CATM{2},MintE2);
    COUP_CATM{3} = cat(2,COUP_CATM{3},MintE3);
    COUP_CATM{4} = cat(2,COUP_CATM{4},MintE4);
end
    
subplot(449)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.5>f/f_p\geq0.8)')
xlabel('Obs E (m^2)')
ylabel('COUP050 E (m^2)')
legend(Bnames,'NumColumns',2)
[~,Sk{9},maxcor{9},rmse{9}] = mod_error(COUP_CATM{1},COUP_CATO{1});
STD{9} = std(COUP_CATM{1});
str = {sprintf('RMSE=%.4f',round(rmse{9}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{9}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{9}.CCall(101),4))};
dim = [0.210560062402496,0.351565030696509,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,10)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.8>f/f_p\geq1.2)')
xlabel('Obs E (m^2)')
ylabel('COUP050 E (m^2)')
[~,Sk{10},maxcor{10},rmse{10}] = mod_error(COUP_CATM{2},COUP_CATO{2});
STD{10} = std(COUP_CATM{2});
str = {sprintf('RMSE=%.4f',round(rmse{10}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{10}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{10}.CCall(101),4))};
dim = [0.408042121684867,0.425949374604002,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,11)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (1.2>f/f_p\geq2.0)')
xlabel('Obs E (m^2)')
ylabel('COUP050 E (m^2)')
[~,Sk{11},maxcor{11},rmse{11}] = mod_error(COUP_CATM{3},COUP_CATO{3});
STD{11} = std(COUP_CATM{3});
str = {sprintf('RMSE=%.4f',round(rmse{11}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{11}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{11}.CCall(101),4))};
dim = [0.609006240249611,0.423420049691664,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,12)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (2.0>f/f_p\geq3.0)')
xlabel('Obs E (m^2)')
ylabel('COUP050 E (m^2)')
[~,Sk{12},maxcor{12},rmse{12}] = mod_error(COUP_CATM{4},COUP_CATO{4});
STD{12} = std(COUP_CATM{4});
str = {sprintf('RMSE=%.4f',round(rmse{12}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{12}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{12}.CCall(101),4))};
dim = [0.804850234009363,0.4222613011401,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN075 Comparison:
    
SWAN075_CATO{1} = []; SWAN075_CATM{1} = [];
SWAN075_CATO{2} = []; SWAN075_CATM{2} = [];
SWAN075_CATO{3} = []; SWAN075_CATM{3} = [];
SWAN075_CATO{4} = []; SWAN075_CATM{4} = [];

for i = 1:9
    eval(sprintf('OintE1 = Buoy.Obs.%s.intE1;',Bnames{i}))
    eval(sprintf('MintE1 = Buoy.SWAN075.%s.intE1;',Bnames{i}))
    eval(sprintf('OintE2 = Buoy.Obs.%s.intE2;',Bnames{i}))
    eval(sprintf('MintE2 = Buoy.SWAN075.%s.intE2;',Bnames{i}))
    eval(sprintf('OintE3 = Buoy.Obs.%s.intE3;',Bnames{i}))
    eval(sprintf('MintE3 = Buoy.SWAN075.%s.intE3;',Bnames{i}))
    eval(sprintf('OintE4 = Buoy.Obs.%s.intE4;',Bnames{i}))
    eval(sprintf('MintE4 = Buoy.SWAN075.%s.intE4;',Bnames{i}))
    
    subplot(4,4,13)
    scatter(OintE1,MintE1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,14)
    scatter(OintE2,MintE2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,15)
    scatter(OintE3,MintE3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,16)
    scatter(OintE4,MintE4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN075_CATO{1} = cat(2,SWAN075_CATO{1},OintE1);
    SWAN075_CATO{2} = cat(2,SWAN075_CATO{2},OintE2);
    SWAN075_CATO{3} = cat(2,SWAN075_CATO{3},OintE3);
    SWAN075_CATO{4} = cat(2,SWAN075_CATO{4},OintE4);
    
    SWAN075_CATM{1} = cat(2,SWAN075_CATM{1},MintE1);
    SWAN075_CATM{2} = cat(2,SWAN075_CATM{2},MintE2);
    SWAN075_CATM{3} = cat(2,SWAN075_CATM{3},MintE3);
    SWAN075_CATM{4} = cat(2,SWAN075_CATM{4},MintE4);
end
    
subplot(4,4,13)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.5>f/f_p\geq0.8)')
xlabel('Obs E (m^2)')
ylabel('SWAN075 E (m^2)')
legend(Bnames,'NumColumns',2)
[~,Sk{13},maxcor{13},rmse{13}] = mod_error(SWAN075_CATM{1},SWAN075_CATO{1});
STD{13} = std(SWAN075_CATM{1});
str = {sprintf('RMSE=%.4f',round(rmse{13}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{13}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{13}.CCall(101),4))};
dim = [0.218560062402496,0.133937174330895,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,14)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.8>f/f_p\geq1.2)')
xlabel('Obs E (m^2)')
ylabel('SWAN075 E (m^2)')
[~,Sk{14},maxcor{14},rmse{14}] = mod_error(SWAN075_CATM{2},SWAN075_CATO{2});
STD{14} = std(SWAN075_CATM{2});
str = {sprintf('RMSE=%.4f',round(rmse{14}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{14}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{14}.CCall(101),4))};
dim = [0.403042121684867,0.203968961111075,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,15)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (1.2>f/f_p\geq2.0)')
xlabel('Obs E (m^2)')
ylabel('SWAN075 E (m^2)')
[~,Sk{15},maxcor{15},rmse{15}] = mod_error(SWAN075_CATM{3},SWAN075_CATO{3});
STD{15} = std(SWAN075_CATM{3});
str = {sprintf('RMSE=%.4f',round(rmse{15}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{15}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{15}.CCall(101),4))};
dim = [0.606006240249611,0.206880332607877,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,16)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (2.0>f/f_p\geq3.0)')
xlabel('Obs E (m^2)')
ylabel('SWAN075 E (m^2)')
[~,Sk{16},maxcor{16},rmse{16}] = mod_error(SWAN075_CATM{4},SWAN075_CATO{4});
STD{16} = std(SWAN075_CATM{4});
str = {sprintf('RMSE=%.4f',round(rmse{16}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{16}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{12}.CCall(101),4))};
dim = [0.807850234009363,0.201369026929001,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % Rename Energy Statistics Variables
E_rmse = rmse;
E_Sk = Sk;
E_maxcor = maxcor;
E_STD = STD;
E_OSTD = OSTD;
clear rmse Sk maxcor STD OSTD


%% Plot Direction Comparison
    
figure(2);clf;
set(gcf,'position',[300,50,1000,940])


    % SWAN050 Comparison:
    
clear SWAN_CATO SWAN_CATM
SWAN_CATO{1} = []; SWAN_CATM{1} = [];
SWAN_CATO{2} = []; SWAN_CATM{2} = [];
SWAN_CATO{3} = []; SWAN_CATM{3} = [];
SWAN_CATO{4} = []; SWAN_CATM{4} = [];    
    
for i = 1:9
    eval(sprintf('O_D_EWM1 = Buoy.Obs.%s.EWMDir1;',Bnames{i}))
    eval(sprintf('M_D_EWM1 = Buoy.SWAN050.%s.EWMDir1;',Bnames{i}))
    eval(sprintf('O_D_EWM2 = Buoy.Obs.%s.EWMDir2;',Bnames{i}))
    eval(sprintf('M_D_EWM2 = Buoy.SWAN050.%s.EWMDir2;',Bnames{i}))
    eval(sprintf('O_D_EWM3 = Buoy.Obs.%s.EWMDir3;',Bnames{i}))
    eval(sprintf('M_D_EWM3 = Buoy.SWAN050.%s.EWMDir3;',Bnames{i}))
    eval(sprintf('O_D_EWM4 = Buoy.Obs.%s.EWMDir4;',Bnames{i}))
    eval(sprintf('M_D_EWM4 = Buoy.SWAN050.%s.EWMDir4;',Bnames{i}))
    
    subplot(441)
    scatter(O_D_EWM1,M_D_EWM1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(442)
    scatter(O_D_EWM2,M_D_EWM2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(443)
    scatter(O_D_EWM3,M_D_EWM3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(444)
    scatter(O_D_EWM4,M_D_EWM4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN_CATO{1} = cat(2,SWAN_CATO{1},O_D_EWM1);
    SWAN_CATO{2} = cat(2,SWAN_CATO{2},O_D_EWM2);
    SWAN_CATO{3} = cat(2,SWAN_CATO{3},O_D_EWM3);
    SWAN_CATO{4} = cat(2,SWAN_CATO{4},O_D_EWM4);
    
    SWAN_CATM{1} = cat(2,SWAN_CATM{1},M_D_EWM1);
    SWAN_CATM{2} = cat(2,SWAN_CATM{2},M_D_EWM2);
    SWAN_CATM{3} = cat(2,SWAN_CATM{3},M_D_EWM3);
    SWAN_CATM{4} = cat(2,SWAN_CATM{4},M_D_EWM4);
end

subplot(441)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.5>f/f_p\geq0.8)')
xlabel('Obs Dir (^o)')
ylabel('SWAN050 Dir (^o)')
legend(Bnames,'fontsize',7,'NumColumns',3,'position',[0.129291915130183,0.872806737902478,0.159126363869018,0.050984934922027])
[~,Sk{1},maxcor{1},rmse{1}] = mod_error(SWAN_CATM{1},SWAN_CATO{1});
STD{1} = std(SWAN_CATM{1});
OSTD{1} = nanstd(SWAN_CATO{1});
str = {sprintf('RMSE=%.4f',round(rmse{1}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{1}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{1}.CCall(101),4))};
dim = [0.211439937597504,0.750050750412623,0.069251170046802,0.057972190034762];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(442)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.8>f/f_p\geq1.2)')
xlabel('Obs Dir (^o)')
ylabel('SWAN050 Dir (^o)')
[~,Sk{2},maxcor{2},rmse{2}] = mod_error(SWAN_CATM{2},SWAN_CATO{2});
STD{2} = std(SWAN_CATM{2});
OSTD{2} = nanstd(SWAN_CATO{2});
str = {sprintf('RMSE=%.4f',round(rmse{2}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{2}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{2}.CCall(101),4))};
dim = [0.424631618361946,0.742669138690294,0.067445363700615,0.063882062344352];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(443)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (1.2>f/f_p\geq2.0)')
xlabel('Obs Dir (^o)')
ylabel('SWAN050 Dir (^o)')
[~,Sk{3},maxcor{3},rmse{3}] = mod_error(SWAN_CATM{3},SWAN_CATO{3});
STD{3} = std(SWAN_CATM{3});
OSTD{3} = nanstd(SWAN_CATO{3});
str = {sprintf('RMSE=%.4f',round(rmse{3}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{3}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{3}.CCall(101),4))};
dim = [0.627549141965681,0.744996310712669,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(444)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (2.0>f/f_p\geq3.0)')
xlabel('Obs Dir (^o)')
ylabel('SWAN050 Dir (^o)')
[~,Sk{4},maxcor{4},rmse{4}] = mod_error(SWAN_CATM{4},SWAN_CATO{4});
STD{4} = nanstd(SWAN_CATM{4});
OSTD{4} = nanstd(SWAN_CATO{4});
str = {sprintf('RMSE=%.4f',round(rmse{4}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{4}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{4}.CCall(101),4))};
dim = [0.836833073322935,0.74622566853397,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN050sm Comparison:

clear SWANsm_CATO SWANsmCATM
SWANsm_CATO{1} = []; SWANsm_CATM{1} = [];
SWANsm_CATO{2} = []; SWANsm_CATM{2} = [];
SWANsm_CATO{3} = []; SWANsm_CATM{3} = [];
SWANsm_CATO{4} = []; SWANsm_CATM{4} = [];

for i = 1:9
    eval(sprintf('O_D_EWM1 = Buoy.Obs.%s.EWMDir1;',Bnames{i}))
    eval(sprintf('M_D_EWM1 = Buoy.SWAN050sm.%s.EWMDir1;',Bnames{i}))
    eval(sprintf('O_D_EWM2 = Buoy.Obs.%s.EWMDir2;',Bnames{i}))
    eval(sprintf('M_D_EWM2 = Buoy.SWAN050sm.%s.EWMDir2;',Bnames{i}))
    eval(sprintf('O_D_EWM3 = Buoy.Obs.%s.EWMDir3;',Bnames{i}))
    eval(sprintf('M_D_EWM3 = Buoy.SWAN050sm.%s.EWMDir3;',Bnames{i}))
    eval(sprintf('O_D_EWM4 = Buoy.Obs.%s.EWMDir4;',Bnames{i}))
    eval(sprintf('M_D_EWM4 = Buoy.SWAN050sm.%s.EWMDir4;',Bnames{i}))
    
    subplot(445)
    scatter(O_D_EWM1,M_D_EWM1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(446)
    scatter(O_D_EWM2,M_D_EWM2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(447)
    scatter(O_D_EWM3,M_D_EWM3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(448)
    scatter(O_D_EWM4,M_D_EWM4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWANsm_CATO{1} = cat(2,SWANsm_CATO{1},O_D_EWM1);
    SWANsm_CATO{2} = cat(2,SWANsm_CATO{2},O_D_EWM2);
    SWANsm_CATO{3} = cat(2,SWANsm_CATO{3},O_D_EWM3);
    SWANsm_CATO{4} = cat(2,SWANsm_CATO{4},O_D_EWM4);
    
    SWANsm_CATM{1} = cat(2,SWANsm_CATM{1},M_D_EWM1);
    SWANsm_CATM{2} = cat(2,SWANsm_CATM{2},M_D_EWM2);
    SWANsm_CATM{3} = cat(2,SWANsm_CATM{3},M_D_EWM3);
    SWANsm_CATM{4} = cat(2,SWANsm_CATM{4},M_D_EWM4);
end

subplot(445)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.5>f/f_p\geq0.8)')
xlabel('Obs Dir (^o)')
ylabel('SWAN050sm Dir (^o)')
legend(Bnames,'NumColumns',3,'fontsize',7,'position',[0.129325065863444,0.654877531276056,0.1615,0.052149987841899])
[~,Sk{5},maxcor{5},rmse{5}] = mod_error(SWANsm_CATM{1},SWANsm_CATO{1});
STD{5} = std(SWANsm_CATM{1});
str = {sprintf('RMSE=%.4f',round(rmse{5}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{5}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{5}.CCall(101),4))};
dim = [0.216439937597504,0.52775621538136,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(446)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.8>f/f_p\geq1.2)')
xlabel('Obs Dir (^o)')
ylabel('SWAN050sm Dir (^o)')
[~,Sk{6},maxcor{6},rmse{6}] = mod_error(SWANsm_CATM{2},SWANsm_CATO{2});
STD{6} = std(SWANsm_CATM{2});
str = {sprintf('RMSE=%.4f',round(rmse{6}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{6}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{6}.CCall(101),4))};
dim = [0.417524180967239,0.526244420481113,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(447)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (1.2>f/f_p\geq2.0)')
xlabel('Obs Dir (^o)')
ylabel('SWAN050sm Dir (^o)')
[~,Sk{7},maxcor{7},rmse{7}] = mod_error(SWANsm_CATM{3},SWANsm_CATO{3});
STD{7} = std(SWANsm_CATM{3});
str = {sprintf('RMSE=%.4f',round(rmse{7}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{7}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{7}.CCall(101),4))};
dim = [0.629748829953199,0.529508838326597,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(448)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (2.0>f/f_p\geq3.0)')
xlabel('Obs Dir (^o)')
ylabel('SWAN050sm Dir (^o)')
[~,Sk{8},maxcor{8},rmse{8}] = mod_error(SWANsm_CATM{4},SWANsm_CATO{4});
STD{8} = std(SWANsm_CATM{4});
str = {sprintf('RMSE=%.4f',round(rmse{8}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{8}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{8}.CCall(101),4))};
dim = [0.830152886115445,0.527261950493204,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % COUP050 Comparison:

clear COUP_CATO COUP_CATM
COUP_CATO{1} = []; COUP_CATM{1} = [];
COUP_CATO{2} = []; COUP_CATM{2} = [];
COUP_CATO{3} = []; COUP_CATM{3} = [];
COUP_CATO{4} = []; COUP_CATM{4} = [];

for i = 1:9
    eval(sprintf('O_D_EWM1 = Buoy.Obs.%s.EWMDir1;',Bnames{i}))
    eval(sprintf('M_D_EWM1 = Buoy.COUP050.%s.EWMDir1;',Bnames{i}))
    eval(sprintf('O_D_EWM2 = Buoy.Obs.%s.EWMDir2;',Bnames{i}))
    eval(sprintf('M_D_EWM2 = Buoy.COUP050.%s.EWMDir2;',Bnames{i}))
    eval(sprintf('O_D_EWM3 = Buoy.Obs.%s.EWMDir3;',Bnames{i}))
    eval(sprintf('M_D_EWM3 = Buoy.COUP050.%s.EWMDir3;',Bnames{i}))
    eval(sprintf('O_D_EWM4 = Buoy.Obs.%s.EWMDir4;',Bnames{i}))
    eval(sprintf('M_D_EWM4 = Buoy.COUP050.%s.EWMDir4;',Bnames{i}))
    
    subplot(449)
    scatter(O_D_EWM1,M_D_EWM1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,10)
    scatter(O_D_EWM2,M_D_EWM2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,11)
    scatter(O_D_EWM3,M_D_EWM3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,12)
    scatter(O_D_EWM4,M_D_EWM4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    COUP_CATO{1} = cat(2,COUP_CATO{1},O_D_EWM1);
    COUP_CATO{2} = cat(2,COUP_CATO{2},O_D_EWM2);
    COUP_CATO{3} = cat(2,COUP_CATO{3},O_D_EWM3);
    COUP_CATO{4} = cat(2,COUP_CATO{4},O_D_EWM4);
    
    COUP_CATM{1} = cat(2,COUP_CATM{1},M_D_EWM1);
    COUP_CATM{2} = cat(2,COUP_CATM{2},M_D_EWM2);
    COUP_CATM{3} = cat(2,COUP_CATM{3},M_D_EWM3);
    COUP_CATM{4} = cat(2,COUP_CATM{4},M_D_EWM4);
end


subplot(449)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.5>f/f_p\geq0.8)')
xlabel('Obs Dir (^o)')
ylabel('COUP050 Dir (^o)')
legend(Bnames,'fontsize',7,'NumColumns',3,'position',[0.131105097064692,0.434936338710707,0.1615,0.052149987841899])
[~,Sk{9},maxcor{9},rmse{9}] = mod_error(COUP_CATM{1},COUP_CATO{1});
STD{9} = std(COUP_CATM{1});
str = {sprintf('RMSE=%.4f',round(rmse{9}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{9}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{9}.CCall(101),4))};
dim = [0.2157,0.309965070729053,0.0698,0.0603];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,10)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.8>f/f_p\geq1.2)')
xlabel('Obs Dir (^o)')
ylabel('COUP050 Dir (^o)')
[~,Sk{10},maxcor{10},rmse{10}] = mod_error(COUP_CATM{2},COUP_CATO{2});
STD{10} = std(COUP_CATM{2});
str = {sprintf('RMSE=%.4f',round(rmse{10}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{10}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{10}.CCall(101),4))};
dim = [0.42486271450858,0.310797003582551,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,11)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (1.2>f/f_p\geq2.0)')
xlabel('Obs Dir (^o)')
ylabel('COUP050 Dir (^o)')
[~,Sk{11},maxcor{11},rmse{11}] = mod_error(COUP_CATM{3},COUP_CATO{3});
STD{11} = std(COUP_CATM{3});
str = {sprintf('RMSE=%.4f',round(rmse{11}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{11}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{11}.CCall(101),4))};
dim = [0.628307332293293,0.310797003582551,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,12)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (2.0>f/f_p\geq3.0)')
xlabel('Obs Dir (^o)')
ylabel('COUP050 Dir (^o)')
[~,Sk{12},maxcor{12},rmse{12}] = mod_error(COUP_CATM{4},COUP_CATO{4});
STD{12} = std(COUP_CATM{4});
str = {sprintf('RMSE=%.4f',round(rmse{12}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{12}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{12}.CCall(101),4))};
dim = [0.833491419656788,0.307391367197594,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN075 Comparison:

clear SWAN075_CATO SWAN075_CATM
SWAN075_CATO{1} = []; SWAN075_CATM{1} = [];
SWAN075_CATO{2} = []; SWAN075_CATM{2} = [];
SWAN075_CATO{3} = []; SWAN075_CATM{3} = [];
SWAN075_CATO{4} = []; SWAN075_CATM{4} = [];

for i = 1:9
    eval(sprintf('O_D_EWM1 = Buoy.Obs.%s.EWMDir1;',Bnames{i}))
    eval(sprintf('M_D_EWM1 = Buoy.SWAN075.%s.EWMDir1;',Bnames{i}))
    eval(sprintf('O_D_EWM2 = Buoy.Obs.%s.EWMDir2;',Bnames{i}))
    eval(sprintf('M_D_EWM2 = Buoy.SWAN075.%s.EWMDir2;',Bnames{i}))
    eval(sprintf('O_D_EWM3 = Buoy.Obs.%s.EWMDir3;',Bnames{i}))
    eval(sprintf('M_D_EWM3 = Buoy.SWAN075.%s.EWMDir3;',Bnames{i}))
    eval(sprintf('O_D_EWM4 = Buoy.Obs.%s.EWMDir4;',Bnames{i}))
    eval(sprintf('M_D_EWM4 = Buoy.SWAN075.%s.EWMDir4;',Bnames{i}))
    
    subplot(4,4,13)
    scatter(O_D_EWM1,M_D_EWM1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,14)
    scatter(O_D_EWM2,M_D_EWM2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,15)
    scatter(O_D_EWM3,M_D_EWM3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,16)
    scatter(O_D_EWM4,M_D_EWM4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN075_CATO{1} = cat(2,SWAN075_CATO{1},O_D_EWM1);
    SWAN075_CATO{2} = cat(2,SWAN075_CATO{2},O_D_EWM2);
    SWAN075_CATO{3} = cat(2,SWAN075_CATO{3},O_D_EWM3);
    SWAN075_CATO{4} = cat(2,SWAN075_CATO{4},O_D_EWM4);
    
    SWAN075_CATM{1} = cat(2,SWAN075_CATM{1},M_D_EWM1);
    SWAN075_CATM{2} = cat(2,SWAN075_CATM{2},M_D_EWM2);
    SWAN075_CATM{3} = cat(2,SWAN075_CATM{3},M_D_EWM3);
    SWAN075_CATM{4} = cat(2,SWAN075_CATM{4},M_D_EWM4);
end


subplot(4,4,13)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.5>f/f_p\geq0.8)')
xlabel('Obs Dir (^o)')
ylabel('SWAN075 Dir (^o)')
legend(Bnames,'fontsize',7,'NumColumns',3,'position',[0.129105097064692,0.219484760908749,0.1615,0.052149987841899])
[~,Sk{13},maxcor{13},rmse{13}] = mod_error(SWAN075_CATM{1},SWAN075_CATO{1});
STD{13} = std(SWAN075_CATM{1});
str = {sprintf('RMSE=%.4f',round(rmse{13}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{13}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{13}.CCall(101),4))};
dim = [0.2167,0.092337214363439,0.0698,0.0603];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,14)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.8>f/f_p\geq1.2)')
xlabel('Obs Dir (^o)')
ylabel('SWAN075 Dir (^o)')
[~,Sk{14},maxcor{14},rmse{14}] = mod_error(SWAN075_CATM{2},SWAN075_CATO{2});
STD{14} = std(SWAN075_CATM{2});
str = {sprintf('RMSE=%.4f',round(rmse{14}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{14}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{14}.CCall(101),4))};
dim = [0.42186271450858,0.089904729371452,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,15)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (1.2>f/f_p\geq2.0)')
xlabel('Obs Dir (^o)')
ylabel('SWAN075 Dir (^o)')
[~,Sk{15},maxcor{15},rmse{15}] = mod_error(SWAN075_CATM{3},SWAN075_CATO{3});
STD{15} = std(SWAN075_CATM{3});
str = {sprintf('RMSE=%.4f',round(rmse{15}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{15}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{15}.CCall(101),4))};
dim = [0.631307332293293,0.093169147216936,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,16)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (2.0>f/f_p\geq3.0)')
xlabel('Obs Dir (^o)')
ylabel('SWAN075 Dir (^o)')
[~,Sk{16},maxcor{16},rmse{16}] = mod_error(SWAN075_CATM{4},SWAN075_CATO{4});
STD{16} = std(SWAN075_CATM{4});
str = {sprintf('RMSE=%.4f',round(rmse{16}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{16}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{16}.CCall(101),4))};
dim = [0.833491419656788,0.090851650113807,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % Rename Direction Statistics Variables
D_rmse = rmse;
D_Sk = Sk;
D_maxcor = maxcor;
D_STD = STD;
D_OSTD = OSTD;
clear rmse Sk maxcor STD OSTD



%% Plot Directional Spread Comparison
    
figure(3);clf; % directional spread figure
set(gcf,'position',[450,50,1000,940])


    % SWAN050 Comparison:
    
clear SWAN_CATO SWAN_CATM
SWAN_CATO{1} = []; SWAN_CATM{1} = [];
SWAN_CATO{2} = []; SWAN_CATM{2} = [];
SWAN_CATO{3} = []; SWAN_CATM{3} = [];
SWAN_CATO{4} = []; SWAN_CATM{4} = [];    
    
for i = 1:9
    eval(sprintf('O_DS_EWM1 = Buoy.Obs.%s.EWMDspread1;',Bnames{i}))
    eval(sprintf('M_DS_EWM1 = Buoy.SWAN050.%s.EWMDspread1;',Bnames{i}))
    eval(sprintf('O_DS_EWM2 = Buoy.Obs.%s.EWMDspread2;',Bnames{i}))
    eval(sprintf('M_DS_EWM2 = Buoy.SWAN050.%s.EWMDspread2;',Bnames{i}))
    eval(sprintf('O_DS_EWM3 = Buoy.Obs.%s.EWMDspread3;',Bnames{i}))
    eval(sprintf('M_DS_EWM3 = Buoy.SWAN050.%s.EWMDspread3;',Bnames{i}))
    eval(sprintf('O_DS_EWM4 = Buoy.Obs.%s.EWMDspread4;',Bnames{i}))
    eval(sprintf('M_DS_EWM4 = Buoy.SWAN050.%s.EWMDspread4;',Bnames{i}))
    
    subplot(441)
    scatter(O_DS_EWM1,M_DS_EWM1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(442)
    scatter(O_DS_EWM2,M_DS_EWM2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(443)
    scatter(O_DS_EWM3,M_DS_EWM3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(444)
    scatter(O_DS_EWM4,M_DS_EWM4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN_CATO{1} = cat(2,SWAN_CATO{1},O_DS_EWM1);
    SWAN_CATO{2} = cat(2,SWAN_CATO{2},O_DS_EWM2);
    SWAN_CATO{3} = cat(2,SWAN_CATO{3},O_DS_EWM3);
    SWAN_CATO{4} = cat(2,SWAN_CATO{4},O_DS_EWM4);
    
    SWAN_CATM{1} = cat(2,SWAN_CATM{1},M_DS_EWM1);
    SWAN_CATM{2} = cat(2,SWAN_CATM{2},M_DS_EWM2);
    SWAN_CATM{3} = cat(2,SWAN_CATM{3},M_DS_EWM3);
    SWAN_CATM{4} = cat(2,SWAN_CATM{4},M_DS_EWM4);
end


subplot(441)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.5>f/f_p\geq0.8)')
xlabel('Obs Spread (^o)')
ylabel('SWAN050 Spread (^o)')
legend(Bnames,'fontsize',7,'NumColumns',3,'location',[0.431105097064692,0.013682318080192,0.170894902935308,0.055958595956804])
[~,Sk{1},maxcor{1},rmse{1}] = mod_error(SWAN_CATM{1},SWAN_CATO{1});
STD{1} = std(SWAN_CATM{1});
OSTD{1} = nanstd(SWAN_CATO{1});
str = {sprintf('RMSE=%.4f',round(rmse{1}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{1}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{1}.CCall(101),4))};
dim = [0.132656786271451,0.867566010210605,0.063790951638066,0.055654692931633];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(442)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (0.8>f/f_p\geq1.2)')
xlabel('Obs Spread (^o)')
ylabel('SWAN050 Spread (^o)')
[~,Sk{2},maxcor{2},rmse{2}] = mod_error(SWAN_CATM{2},SWAN_CATO{2});
STD{2} = std(SWAN_CATM{2});
OSTD{2} = nanstd(SWAN_CATO{2});
str = {sprintf('RMSE=%.4f',round(rmse{2}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{2}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{2}.CCall(101),4))};
dim = [0.336408529438389,0.859785960466193,0.067445363700615,0.063882062344352];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(443)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (1.2>f/f_p\geq2.0)')
xlabel('Obs Spread (^o)')
ylabel('SWAN050 Spread (^o)')
[~,Sk{3},maxcor{3},rmse{3}] = mod_error(SWAN_CATM{3},SWAN_CATO{3});
STD{3} = std(SWAN_CATM{3});
OSTD{3} = nanstd(SWAN_CATO{3});
str = {sprintf('RMSE=%.4f',round(rmse{3}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{3}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{3}.CCall(101),4))};
dim = [0.628549141965681,0.866889345234297,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(444)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050 (2.0>f/f_p\geq3.0)')
xlabel('Obs Spread (^o)')
ylabel('SWAN050 Spread (^o)')
[~,Sk{4},maxcor{4},rmse{4}] = mod_error(SWAN_CATM{4},SWAN_CATO{4});
STD{4} = nanstd(SWAN_CATM{4});
OSTD{4} = nanstd(SWAN_CATO{4});
str = {sprintf('RMSE=%.4f',round(rmse{4}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{4}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{4}.CCall(101),4))};
dim = [0.840093603744152,0.866889345234298,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN050sm Comparison:

clear SWANsm_CATO SWANsmCATM
SWANsm_CATO{1} = []; SWANsm_CATM{1} = [];
SWANsm_CATO{2} = []; SWANsm_CATM{2} = [];
SWANsm_CATO{3} = []; SWANsm_CATM{3} = [];
SWANsm_CATO{4} = []; SWANsm_CATM{4} = [];

for i = 1:9
    eval(sprintf('O_DS_EWM1 = Buoy.Obs.%s.EWMDspread1;',Bnames{i}))
    eval(sprintf('M_DS_EWM1 = Buoy.SWAN050sm.%s.EWMDspread1;',Bnames{i}))
    eval(sprintf('O_DS_EWM2 = Buoy.Obs.%s.EWMDspread2;',Bnames{i}))
    eval(sprintf('M_DS_EWM2 = Buoy.SWAN050sm.%s.EWMDspread2;',Bnames{i}))
    eval(sprintf('O_DS_EWM3 = Buoy.Obs.%s.EWMDspread3;',Bnames{i}))
    eval(sprintf('M_DS_EWM3 = Buoy.SWAN050sm.%s.EWMDspread3;',Bnames{i}))
    eval(sprintf('O_DS_EWM4 = Buoy.Obs.%s.EWMDspread4;',Bnames{i}))
    eval(sprintf('M_DS_EWM4 = Buoy.SWAN050sm.%s.EWMDspread4;',Bnames{i}))
    
    subplot(445)
    scatter(O_DS_EWM1,M_DS_EWM1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(446)
    scatter(O_DS_EWM2,M_DS_EWM2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(447)
    scatter(O_DS_EWM3,M_DS_EWM3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(448)
    scatter(O_DS_EWM4,M_DS_EWM4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWANsm_CATO{1} = cat(2,SWANsm_CATO{1},O_DS_EWM1);
    SWANsm_CATO{2} = cat(2,SWANsm_CATO{2},O_DS_EWM2);
    SWANsm_CATO{3} = cat(2,SWANsm_CATO{3},O_DS_EWM3);
    SWANsm_CATO{4} = cat(2,SWANsm_CATO{4},O_DS_EWM4);
    
    SWANsm_CATM{1} = cat(2,SWANsm_CATM{1},M_DS_EWM1);
    SWANsm_CATM{2} = cat(2,SWANsm_CATM{2},M_DS_EWM2);
    SWANsm_CATM{3} = cat(2,SWANsm_CATM{3},M_DS_EWM3);
    SWANsm_CATM{4} = cat(2,SWANsm_CATM{4},M_DS_EWM4);
end

subplot(445)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.5>f/f_p\geq0.8)')
xlabel('Obs Spread (^o)')
ylabel('SWAN050sm Spread (^o)')
%legend(Bnames,'fontsize',7,'NumColumns',3,'location',[0.128511883928935,0.573849611598886,0.159126363869018,0.050984934922027])
[~,Sk{5},maxcor{5},rmse{5}] = mod_error(SWANsm_CATM{1},SWANsm_CATO{1});
STD{5} = std(SWANsm_CATM{1});
str = {sprintf('RMSE=%.4f',round(rmse{5}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{5}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{5}.CCall(101),4))};
dim = [0.130776911076443,0.642289494412802,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(446)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (0.8>f/f_p\geq1.2)')
xlabel('Obs Spread (^o)')
ylabel('SWAN050sm Spread (^o)')
[~,Sk{6},maxcor{6},rmse{6}] = mod_error(SWANsm_CATM{2},SWANsm_CATO{2});
STD{6} = std(SWANsm_CATM{2});
str = {sprintf('RMSE=%.4f',round(rmse{6}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{6}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{6}.CCall(101),4))};
dim = [0.341421216848675,0.640831917344675,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(447)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (1.2>f/f_p\geq2.0)')
xlabel('Obs Spread (^o)')
ylabel('SWAN050sm Spread (^o)')
[~,Sk{7},maxcor{7},rmse{7}] = mod_error(SWANsm_CATM{3},SWANsm_CATO{3});
STD{7} = std(SWANsm_CATM{3});
str = {sprintf('RMSE=%.4f',round(rmse{7}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{7}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{7}.CCall(101),4))};
dim = [0.546945397815913,0.643502460796484,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(448)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN050sm (2.0>f/f_p\geq3.0)')
xlabel('Obs Spread (^o)')
ylabel('SWAN050sm Spread (^o)')
[~,Sk{8},maxcor{8},rmse{8}] = mod_error(SWANsm_CATM{4},SWANsm_CATO{4});
STD{8} = std(SWANsm_CATM{4});
str = {sprintf('RMSE=%.4f',round(rmse{8}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{8}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{8}.CCall(101),4))};
dim = [0.751909516380656,0.64567873936014,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % COUP050 Comparison:

clear COUP_CATO COUP_CATM
COUP_CATO{1} = []; COUP_CATM{1} = [];
COUP_CATO{2} = []; COUP_CATM{2} = [];
COUP_CATO{3} = []; COUP_CATM{3} = [];
COUP_CATO{4} = []; COUP_CATM{4} = [];

for i = 1:9
    eval(sprintf('O_DS_EWM1 = Buoy.Obs.%s.EWMDspread1;',Bnames{i}))
    eval(sprintf('M_DS_EWM1 = Buoy.COUP050.%s.EWMDspread1;',Bnames{i}))
    eval(sprintf('O_DS_EWM2 = Buoy.Obs.%s.EWMDspread2;',Bnames{i}))
    eval(sprintf('M_DS_EWM2 = Buoy.COUP050.%s.EWMDspread2;',Bnames{i}))
    eval(sprintf('O_DS_EWM3 = Buoy.Obs.%s.EWMDspread3;',Bnames{i}))
    eval(sprintf('M_DS_EWM3 = Buoy.COUP050.%s.EWMDspread3;',Bnames{i}))
    eval(sprintf('O_DS_EWM4 = Buoy.Obs.%s.EWMDspread4;',Bnames{i}))
    eval(sprintf('M_DS_EWM4 = Buoy.COUP050.%s.EWMDspread4;',Bnames{i}))
    
    subplot(449)
    scatter(O_DS_EWM1,M_DS_EWM1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,10)
    scatter(O_DS_EWM2,M_DS_EWM2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,11)
    scatter(O_DS_EWM3,M_DS_EWM3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,12)
    scatter(O_DS_EWM4,M_DS_EWM4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    COUP_CATO{1} = cat(2,COUP_CATO{1},O_DS_EWM1);
    COUP_CATO{2} = cat(2,COUP_CATO{2},O_DS_EWM2);
    COUP_CATO{3} = cat(2,COUP_CATO{3},O_DS_EWM3);
    COUP_CATO{4} = cat(2,COUP_CATO{4},O_DS_EWM4);
    
    COUP_CATM{1} = cat(2,COUP_CATM{1},M_DS_EWM1);
    COUP_CATM{2} = cat(2,COUP_CATM{2},M_DS_EWM2);
    COUP_CATM{3} = cat(2,COUP_CATM{3},M_DS_EWM3);
    COUP_CATM{4} = cat(2,COUP_CATM{4},M_DS_EWM4);
end

subplot(449)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.5>f/f_p\geq0.8)')
xlabel('Obs Spread (^o)')
ylabel('COUP050 Spread (^o)')
%legend(Bnames,'fontsize',7,'NumColumns',3,'location',[0.129291915130183,0.27338611217826,0.159126363869018,0.050984934922027])
[~,Sk{9},maxcor{9},rmse{9}] = mod_error(COUP_CATM{1},COUP_CATO{1});
STD{9} = std(COUP_CATM{1});
str = {sprintf('RMSE=%.4f',round(rmse{9}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{9}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{9}.CCall(101),4))};
dim = [0.134876755070203,0.420242632553533,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,10)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (0.8>f/f_p\geq1.2)')
xlabel('Obs Spread (^o)')
ylabel('COUP050 Spread (^o)')
[~,Sk{10},maxcor{10},rmse{10}] = mod_error(COUP_CATM{2},COUP_CATO{2});
STD{10} = std(COUP_CATM{2});
str = {sprintf('RMSE=%.4f',round(rmse{10}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{10}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{10}.CCall(101),4))};
dim = [0.339639625585023,0.420296850385653,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,11)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (1.2>f/f_p\geq2.0)')
xlabel('Obs Spread (^o)')
ylabel('COUP050 Spread (^o)')
[~,Sk{11},maxcor{11},rmse{11}] = mod_error(COUP_CATM{3},COUP_CATO{3});
STD{11} = std(COUP_CATM{3});
str = {sprintf('RMSE=%.4f',round(rmse{11}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{11}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{11}.CCall(101),4))};
dim = [0.545703588143527,0.425808156064529,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,12)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs COUP050 (2.0>f/f_p\geq3.0)')
xlabel('Obs Spread (^o)')
ylabel('COUP050 Spread (^o)')
[~,Sk{12},maxcor{12},rmse{12}] = mod_error(COUP_CATM{4},COUP_CATO{4});
STD{12} = std(COUP_CATM{4});
str = {sprintf('RMSE=%.4f',round(rmse{12}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{12}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{12}.CCall(101),4))};
dim = [0.751687987519502,0.422473128949309,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % SWAN075 Comparison:

clear SWAN075_CATO SWAN075_CATM
SWAN075_CATO{1} = []; SWAN075_CATM{1} = [];
SWAN075_CATO{2} = []; SWAN075_CATM{2} = [];
SWAN075_CATO{3} = []; SWAN075_CATM{3} = [];
SWAN075_CATO{4} = []; SWAN075_CATM{4} = [];

for i = 1:9
    eval(sprintf('O_DS_EWM1 = Buoy.Obs.%s.EWMDspread1;',Bnames{i}))
    eval(sprintf('M_DS_EWM1 = Buoy.SWAN075.%s.EWMDspread1;',Bnames{i}))
    eval(sprintf('O_DS_EWM2 = Buoy.Obs.%s.EWMDspread2;',Bnames{i}))
    eval(sprintf('M_DS_EWM2 = Buoy.SWAN075.%s.EWMDspread2;',Bnames{i}))
    eval(sprintf('O_DS_EWM3 = Buoy.Obs.%s.EWMDspread3;',Bnames{i}))
    eval(sprintf('M_DS_EWM3 = Buoy.SWAN075.%s.EWMDspread3;',Bnames{i}))
    eval(sprintf('O_DS_EWM4 = Buoy.Obs.%s.EWMDspread4;',Bnames{i}))
    eval(sprintf('M_DS_EWM4 = Buoy.SWAN075.%s.EWMDspread4;',Bnames{i}))
    
    subplot(4,4,13)
    scatter(O_DS_EWM1,M_DS_EWM1,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,14)
    scatter(O_DS_EWM2,M_DS_EWM2,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,15)
    scatter(O_DS_EWM3,M_DS_EWM3,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    subplot(4,4,16)
    scatter(O_DS_EWM4,M_DS_EWM4,MSize,pformat{i},'LineWidth',1.5,'MarkerEdgeAlpha',transformat(i))
    hold on
    
    SWAN075_CATO{1} = cat(2,SWAN075_CATO{1},O_DS_EWM1);
    SWAN075_CATO{2} = cat(2,SWAN075_CATO{2},O_DS_EWM2);
    SWAN075_CATO{3} = cat(2,SWAN075_CATO{3},O_DS_EWM3);
    SWAN075_CATO{4} = cat(2,SWAN075_CATO{4},O_DS_EWM4);
    
    SWAN075_CATM{1} = cat(2,SWAN075_CATM{1},M_DS_EWM1);
    SWAN075_CATM{2} = cat(2,SWAN075_CATM{2},M_DS_EWM2);
    SWAN075_CATM{3} = cat(2,SWAN075_CATM{3},M_DS_EWM3);
    SWAN075_CATM{4} = cat(2,SWAN075_CATM{4},M_DS_EWM4);
end

subplot(4,4,13)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.5>f/f_p\geq0.8)')
xlabel('Obs Spread (^o)')
ylabel('SWAN075 Spread (^o)')
%legend(Bnames,'fontsize',7,'NumColumns',3,'location',[0.129291915130183,0.27338611217826,0.159126363869018,0.050984934922027])
[~,Sk{13},maxcor{13},rmse{13}] = mod_error(SWAN075_CATM{1},SWAN075_CATO{1});
STD{13} = std(SWAN075_CATM{1});
str = {sprintf('RMSE=%.4f',round(rmse{13}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{13}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{13}.CCall(101),4))};
dim = [0.130876755070203,0.203702915469747,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,14)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (0.8>f/f_p\geq1.2)')
xlabel('Obs Spread (^o)')
ylabel('SWAN075 Spread (^o)')
[~,Sk{14},maxcor{14},rmse{14}] = mod_error(SWAN075_CATM{2},SWAN075_CATO{2});
STD{14} = std(SWAN075_CATM{2});
str = {sprintf('RMSE=%.4f',round(rmse{14}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{14}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{14}.CCall(101),4))};
dim = [0.338639625585023,0.203757133301866,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,15)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (1.2>f/f_p\geq2.0)')
xlabel('Obs Spread (^o)')
ylabel('SWAN075 Spread (^o)')
[~,Sk{15},maxcor{15},rmse{15}] = mod_error(SWAN075_CATM{3},SWAN075_CATO{3});
STD{15} = std(SWAN075_CATM{3});
str = {sprintf('RMSE=%.4f',round(rmse{15}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{15}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{15}.CCall(101),4))};
dim = [0.545703588143527,0.206004021135258,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')

subplot(4,4,16)
a = gca; xL = a.XLim; yL = a.YLim;
set(a,'color',lightgrey)
plot(one_to_one,one_to_one,'--k')
xlim([xL]); ylim([yL])
axis equal; grid on;
title('Obs vs SWAN075 (2.0>f/f_p\geq3.0)')
xlabel('Obs Spread (^o)')
ylabel('SWAN075 Spread (^o)')
[~,Sk{16},maxcor{16},rmse{16}] = mod_error(SWAN075_CATM{4},SWAN075_CATO{4});
STD{16} = std(SWAN075_CATM{4});
str = {sprintf('RMSE=%.4f',round(rmse{16}.r1,4)),...
    sprintf('Sk=%.4f',round(Sk{16}.sk1,4)),...
    sprintf('CC=%.4f',round(maxcor{16}.CCall(101),4))};
dim = [0.749687987519502,0.204845272583694,0.069812790663585,0.060254923230941];
annotation('textbox',dim,'String',str,'fontsize',8,'EdgeColor','none')


    % Rename Directional Spread Statistics Variables
DS_rmse = rmse;
DS_Sk = Sk;
DS_maxcor = maxcor;
DS_STD = STD;
DS_OSTD = OSTD;
clear rmse Sk maxcor STD OSTD


%% Calculate Statistics for Obs and Models for each of the Buoys

% Preallocate cells (to make all same size for addtion):
Ei_Sk1 = zeros(16,9);Ei_Sk2 = zeros(16,9);Ei_Sk3 = zeros(16,9);Ei_Sk4 = zeros(16,9);
Ei_maxcor1 = zeros(16,9);Ei_maxcor2 = zeros(16,9);Ei_maxcor3 = zeros(16,9);Ei_maxcor4 = zeros(16,9);
Ei_rmse1 = zeros(16,9);Ei_rmse2 = zeros(16,9);Ei_rmse3 = zeros(16,9);Ei_rmse4 = zeros(16,9);
Eistd_M1 = zeros(16,9);Eistd_M2 = zeros(16,9);Eistd_M3 = zeros(16,9);Eistd_M4 = zeros(16,9);
Eistd_O1 = zeros(16,9);Eistd_O2 = zeros(16,9);Eistd_O3 = zeros(16,9);Eistd_O4 = zeros(16,9);
Di_Sk1 = zeros(16,9);Di_Sk2 = zeros(16,9);Di_Sk3 = zeros(16,9);Di_Sk4 = zeros(16,9);
Di_maxcor1 = zeros(16,9);Di_maxcor2 = zeros(16,9);Di_maxcor3 = zeros(16,9);Di_maxcor4 = zeros(16,9);
Di_rmse1 = zeros(16,9);Di_rmse2 = zeros(16,9);Di_rmse3 = zeros(16,9);Di_rmse4 = zeros(16,9);
Distd_M1 = zeros(16,9);Distd_M2 = zeros(16,9);Distd_M3 = zeros(16,9);Distd_M4 = zeros(16,9);
Distd_O1 = zeros(16,9);Distd_O2 = zeros(16,9);Distd_O3 = zeros(16,9);Distd_O4 = zeros(16,9);
DSi_Sk1 = zeros(16,9);DSi_Sk2 = zeros(16,9);DSi_Sk3 = zeros(16,9);DSi_Sk4 = zeros(16,9);
DSi_maxcor1 = zeros(16,9);DSi_maxcor2 = zeros(16,9);DSi_maxcor3 = zeros(16,9);DSi_maxcor4 = zeros(16,9);
DSi_rmse1 = zeros(16,9);DSi_rmse2 = zeros(16,9);DSi_rmse3 = zeros(16,9);DSi_rmse4 = zeros(16,9);
DSistd_M1 = zeros(16,9);DSistd_M2 = zeros(16,9);DSistd_M3 = zeros(16,9);DSistd_M4 = zeros(16,9);
DSistd_O1 = zeros(16,9);DSistd_O2 = zeros(16,9);DSistd_O3 = zeros(16,9);DSistd_O4 = zeros(16,9);


for MO = 1:4 % for the 4 models
    
    for j = 1:9 % for the 9 buoys
        
        if MO == 1
            eval(sprintf('Mod = Buoy.SWAN050.%s;',Bnames{j})) % SWAN050
        elseif MO == 2
            eval(sprintf('Mod = Buoy.SWAN050sm.%s;',Bnames{j})) % SWAN050sm
        elseif MO == 3
            eval(sprintf('Mod = Buoy.COUP050.%s;',Bnames{j})) % COUP050
        else
            eval(sprintf('Mod = Buoy.SWAN075.%s;',Bnames{j})) % SWAN075
        end
        
        eval(sprintf('Obs = Buoy.Obs.%s;',Bnames{j})) % observed
        
        % Energy (4 freq bands)
        [~,Sk,maxcor,rmse] = mod_error(Mod.intE1,Obs.intE1);
        Ei_Sk1(MO*4-3,j) = Sk.sk1; 
        Ei_maxcor1(MO*4-3,j) = maxcor.CCall(101); 
        Ei_rmse1(MO*4-3,j) = rmse.r1;
        Eistd_M1(MO*4-3,j) = nanstd(Mod.intE1);
        Eistd_O1(MO*4-3,j) = nanstd(Obs.intE1);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.intE2,Obs.intE2);
        Ei_Sk2(MO*4-2,j) = Sk.sk1; 
        Ei_maxcor2(MO*4-2,j) = maxcor.CCall(101); 
        Ei_rmse2(MO*4-2,j) = rmse.r1;
        Eistd_M2(MO*4-2,j) = nanstd(Mod.intE2);
        Eistd_O2(MO*4-2,j) = nanstd(Obs.intE2);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.intE3,Obs.intE3);
        Ei_Sk3(MO*4-1,j) = Sk.sk1; 
        Ei_maxcor3(MO*4-1,j) = maxcor.CCall(101); 
        Ei_rmse3(MO*4-1,j) = rmse.r1;
        Eistd_M3(MO*4-1,j) = nanstd(Mod.intE3);
        Eistd_O3(MO*4-1,j) = nanstd(Obs.intE3);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.intE4,Obs.intE4);
        Ei_Sk4(MO*4,j) = Sk.sk1; 
        Ei_maxcor4(MO*4,j) = maxcor.CCall(101); 
        Ei_rmse4(MO*4,j) = rmse.r1;
        Eistd_M4(MO*4,j) = nanstd(Mod.intE4);
        Eistd_O4(MO*4,j) = nanstd(Obs.intE4);
        
        
        % Direction (4 freq bands)
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMDir1,Obs.EWMDir1);
        Di_Sk1(MO*4-3,j) = Sk.sk1; 
        Di_maxcor1(MO*4-3,j) = maxcor.CCall(101); 
        Di_rmse1(MO*4-3,j) = rmse.r1;
        Distd_M1(MO*4-3,j) = nanstd(Mod.EWMDir1);
        Distd_O1(MO*4-3,j) = nanstd(Obs.EWMDir1);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMDir2,Obs.EWMDir2);
        Di_Sk2(MO*4-2,j) = Sk.sk1; 
        Di_maxcor2(MO*4-2,j) = maxcor.CCall(101); 
        Di_rmse2(MO*4-2,j) = rmse.r1;
        Distd_M2(MO*4-2,j) = nanstd(Mod.EWMDir2);
        Distd_O2(MO*4-2,j) = nanstd(Obs.EWMDir2);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMDir3,Obs.EWMDir3);
        Di_Sk3(MO*4-1,j) = Sk.sk1; 
        Di_maxcor3(MO*4-1,j) = maxcor.CCall(101); 
        Di_rmse3(MO*4-1,j) = rmse.r1;
        Distd_M3(MO*4-1,j) = nanstd(Mod.EWMDir3);
        Distd_O3(MO*4-1,j) = nanstd(Obs.EWMDir3);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMDir4,Obs.EWMDir4);
        Di_Sk4(MO*4,j) = Sk.sk1; 
        Di_maxcor4(MO*4,j) = maxcor.CCall(101); 
        Di_rmse4(MO*4,j) = rmse.r1;
        Distd_M4(MO*4,j) = nanstd(Mod.EWMDir4);
        Distd_O4(MO*4,j) = nanstd(Obs.EWMDir4);
            
        
        % Directional Spread (4 freq bands)
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMDspread1,Obs.EWMDspread1);
        DSi_Sk1(MO*4-3,j) = Sk.sk1; 
        DSi_maxcor1(MO*4-3,j) = maxcor.CCall(101); 
        DSi_rmse1(MO*4-3,j) = rmse.r1;
        DSistd_M1(MO*4-3,j) = nanstd(Mod.EWMDspread1);
        DSistd_O1(MO*4-3,j) = nanstd(Obs.EWMDspread1);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMDspread2,Obs.EWMDspread2);
        DSi_Sk2(MO*4-2,j) = Sk.sk1; 
        DSi_maxcor2(MO*4-2,j) = maxcor.CCall(101); 
        DSi_rmse2(MO*4-2,j) = rmse.r1;
        DSistd_M2(MO*4-2,j) = nanstd(Mod.EWMDspread2);
        DSistd_O2(MO*4-2,j) = nanstd(Obs.EWMDspread2);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMDspread3,Obs.EWMDspread3);
        DSi_Sk3(MO*4-1,j) = Sk.sk1; 
        DSi_maxcor3(MO*4-1,j) = maxcor.CCall(101); 
        DSi_rmse3(MO*4-1,j) = rmse.r1;
        DSistd_M3(MO*4-1,j) = nanstd(Mod.EWMDspread3);
        DSistd_O3(MO*4-1,j) = nanstd(Obs.EWMDspread3);
        
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMDspread4,Obs.EWMDspread4);
        DSi_Sk4(MO*4,j) = Sk.sk1; 
        DSi_maxcor4(MO*4,j) = maxcor.CCall(101); 
        DSi_rmse4(MO*4,j) = rmse.r1;
        DSistd_M4(MO*4,j) = nanstd(Mod.EWMDspread4);
        DSistd_O4(MO*4,j) = nanstd(Obs.EWMDspread4);
        
    end
end

% Add all together (to create one matrix with all)
Ei_Sk = Ei_Sk1 + Ei_Sk2 + Ei_Sk3 + Ei_Sk4;
Ei_maxcor = Ei_maxcor1 + Ei_maxcor2 + Ei_maxcor3 + Ei_maxcor4;
Ei_rmse = Ei_rmse1 + Ei_rmse2 + Ei_rmse3 + Ei_rmse4;
Ei_std_M = Eistd_M1 + Eistd_M2 + Eistd_M3 + Eistd_M4;
Ei_std_O = Eistd_O1 + Eistd_O2 + Eistd_O3 + Eistd_O4;
Di_Sk = Di_Sk1 + Di_Sk2 + Di_Sk3 + Di_Sk4;
Di_maxcor = Di_maxcor1 + Di_maxcor2 + Di_maxcor3 + Di_maxcor4;
Di_rmse = Di_rmse1 + Di_rmse2 + Di_rmse3 + Di_rmse4;
Di_std_M = Distd_M1 + Distd_M2 + Distd_M3 + Distd_M4;
Di_std_O = Distd_O1 + Distd_O2 + Distd_O3 + Distd_O4;
DSi_Sk = DSi_Sk1 + DSi_Sk2 + DSi_Sk3 + DSi_Sk4;
DSi_maxcor = DSi_maxcor1 + DSi_maxcor2 + DSi_maxcor3 + DSi_maxcor4;
DSi_rmse = DSi_rmse1 + DSi_rmse2 + DSi_rmse3 + DSi_rmse4;
DSi_std_M = DSistd_M1 + DSistd_M2 + DSistd_M3 + DSistd_M4;
DSi_std_O = DSistd_O1 + DSistd_O2 + DSistd_O3 + DSistd_O4;

% Clear up individual frequency variables
clear Ei_Sk1 Ei_Sk2 Ei_Sk3 Ei_Sk3 Ei_maxcor1 + Ei_maxcor2 + Ei_maxcor3 + Ei_maxcor4
clear Ei_rmse1 + Ei_rmse2 + Ei_rmse3 + Ei_rmse4 Eistd_M1 + Eistd_M2 + Eistd_M3 + Eistd_M4
clear Eistd_O1 + Eistd_O2 + Eistd_O3 + Eistd_O4 Di_Sk1 + Di_Sk2 + Di_Sk3 + Di_Sk4
clear Di_maxcor1 + Di_maxcor2 + Di_maxcor3 + Di_maxcor4 Di_rmse1 + Di_rmse2 + Di_rmse3 + Di_rmse4
clear Distd_M1 + Distd_M2 + Distd_M3 + Distd_M4 Distd_O1 + Distd_O2 + Distd_O3 + Distd_O4
clear DSi_Sk1 + DSi_Sk2 + DSi_Sk3 + DSi_Sk4 DSi_maxcor1 + DSi_maxcor2 + DSi_maxcor3 + DSi_maxcor4
clear DSi_rmse1 + DSi_rmse2 + DSi_rmse3 + DSi_rmse4 DSistd_M1 + DSistd_M2 + DSistd_M3 + DSistd_M4
clear DSistd_O1 + DSistd_O2 + DSistd_O3 + DSistd_O4 j MO


%% Statistics Tables

        
Model = ["SWAN050";"SWAN050";"SWAN050";"SWAN050";"SWAN050sm";"SWAN050sm";...
    "SWAN050sm";"SWAN050sm";"COUP050";"COUP050";"COUP050";"COUP050"...
    ;"SWAN075";"SWAN075";"SWAN075";"SWAN075"];
f_Band = ["f/fp=0.5-0.8";"f/fp=0.8-1.2";"f/fp=1.2-2.0";"f/fp=2.0-3.0";...
    "f/fp=0.5-0.8";"f/fp=0.8-1.2";"f/fp=1.2-2.0";"f/fp=2.0-3.0";...
    "f/fp=0.5-0.8";"f/fp=0.8-1.2";"f/fp=1.2-2.0";"f/fp=2.0-3.0";...
    "f/fp=0.5-0.8";"f/fp=0.8-1.2";"f/fp=1.2-2.0";"f/fp=2.0-3.0"];


    % ENERGY TABLE:
RMSE = [E_rmse{1}.r1;E_rmse{2}.r1;E_rmse{3}.r1;E_rmse{4}.r1;E_rmse{5}.r1;...
    E_rmse{6}.r1;E_rmse{7}.r1;E_rmse{8}.r1;E_rmse{9}.r1;E_rmse{10}.r1;...
    E_rmse{11}.r1;E_rmse{12}.r1;E_rmse{13}.r1;E_rmse{14}.r1;...
    E_rmse{15}.r1;E_rmse{16}.r1];
RMSE = cat(2,RMSE,Ei_rmse);
RMSE = round(RMSE,3);

SS = [E_Sk{1}.sk1;E_Sk{2}.sk1;E_Sk{3}.sk1;E_Sk{4}.sk1;E_Sk{5}.sk1;E_Sk{6}.sk1;...
    E_Sk{7}.sk1;E_Sk{8}.sk1;E_Sk{9}.sk1;E_Sk{10}.sk1;E_Sk{11}.sk1;...
    E_Sk{12}.sk1;E_Sk{13}.sk1;E_Sk{14}.sk1;E_Sk{15}.sk1;E_Sk{16}.sk1];
SS = cat(2,SS,Ei_Sk);
SS = round(SS,3);

MaxCor = [E_maxcor{1}.CCall(101);E_maxcor{2}.CCall(101);E_maxcor{3}.CCall(101);...
    E_maxcor{4}.CCall(101);E_maxcor{5}.CCall(101);E_maxcor{6}.CCall(101);...
    E_maxcor{7}.CCall(101);E_maxcor{8}.CCall(101);E_maxcor{9}.CCall(101);...
    E_maxcor{10}.CCall(101);E_maxcor{11}.CCall(101);E_maxcor{12}.CCall(101)...
    ;E_maxcor{13}.CCall(101);E_maxcor{14}.CCall(101);E_maxcor{15}.CCall(101);...
    E_maxcor{16}.CCall(101)];
MaxCor = cat(2,MaxCor,Ei_maxcor);
MaxCor = round(MaxCor,3);

STD = [E_STD{1};E_STD{2};E_STD{3};E_STD{4};E_STD{5};E_STD{6};E_STD{7};...
    E_STD{8};E_STD{9};E_STD{10};E_STD{11};E_STD{12};E_STD{13};E_STD{14}...
    ;E_STD{15};E_STD{16}];
STD = cat(2,STD,Ei_std_M);
STD = round(STD,3);

EnergyStats = table(Model,f_Band,RMSE,SS,MaxCor,STD)

clear RMSE SS MaxCor STD


    % DIRECTION TABLE:
RMSE = [D_rmse{1}.r1;D_rmse{2}.r1;D_rmse{3}.r1;D_rmse{4}.r1;D_rmse{5}.r1;...
    D_rmse{6}.r1;D_rmse{7}.r1;D_rmse{8}.r1;D_rmse{9}.r1;D_rmse{10}.r1;...
    D_rmse{11}.r1;D_rmse{12}.r1;D_rmse{13}.r1;D_rmse{14}.r1;D_rmse{15}.r1...
    ;D_rmse{16}.r1];
RMSE = cat(2,RMSE,Di_rmse);
RMSE = round(RMSE,3);

SS = [D_Sk{1}.sk1;D_Sk{2}.sk1;D_Sk{3}.sk1;D_Sk{4}.sk1;D_Sk{5}.sk1;D_Sk{6}.sk1;...
    D_Sk{7}.sk1;D_Sk{8}.sk1;D_Sk{9}.sk1;D_Sk{10}.sk1;D_Sk{11}.sk1;D_Sk{12}.sk1...
    ;D_Sk{13}.sk1;D_Sk{14}.sk1;D_Sk{15}.sk1;D_Sk{16}.sk1];
SS = cat(2,SS,Di_Sk);
SS = round(SS,3);

MaxCor = [D_maxcor{1}.CCall(101);D_maxcor{2}.CCall(101);D_maxcor{3}.CCall(101);...
    D_maxcor{4}.CCall(101);D_maxcor{5}.CCall(101);D_maxcor{6}.CCall(101);...
    D_maxcor{7}.CCall(101);D_maxcor{8}.CCall(101);D_maxcor{9}.CCall(101);...
    D_maxcor{10}.CCall(101);D_maxcor{11}.CCall(101);D_maxcor{12}.CCall(101)...
    ;D_maxcor{13}.CCall(101);D_maxcor{14}.CCall(101);D_maxcor{15}.CCall(101)...
    ;D_maxcor{16}.CCall(101)];
MaxCor = cat(2,MaxCor,Di_maxcor);
MaxCor = round(MaxCor,3);

STD = [D_STD{1};D_STD{2};D_STD{3};D_STD{4};D_STD{5};D_STD{6};D_STD{7};...
    D_STD{8};D_STD{9};D_STD{10};D_STD{11};D_STD{12};D_STD{13};D_STD{14}...
    ;D_STD{15};D_STD{16}];
STD = cat(2,STD,Di_std_M);
STD = round(STD,3);

DirectionStats = table(Model,f_Band,RMSE,SS,MaxCor,STD)

clear RMSE SS MaxCor STD


% DIRECTIONAL SPREAD TABLE:
    % DIRECTION TABLE:
RMSE = [DS_rmse{1}.r1;DS_rmse{2}.r1;DS_rmse{3}.r1;DS_rmse{4}.r1;DS_rmse{5}.r1;...
    DS_rmse{6}.r1;DS_rmse{7}.r1;DS_rmse{8}.r1;DS_rmse{9}.r1;DS_rmse{10}.r1;...
    DS_rmse{11}.r1;DS_rmse{12}.r1;DS_rmse{13}.r1;DS_rmse{14}.r1;DS_rmse{15}.r1...
    ;DS_rmse{16}.r1];
RMSE = cat(2,RMSE,DSi_rmse);
RMSE = round(RMSE,3);

SS = [DS_Sk{1}.sk1;DS_Sk{2}.sk1;DS_Sk{3}.sk1;DS_Sk{4}.sk1;DS_Sk{5}.sk1;DS_Sk{6}.sk1;...
    DS_Sk{7}.sk1;DS_Sk{8}.sk1;DS_Sk{9}.sk1;DS_Sk{10}.sk1;DS_Sk{11}.sk1;DS_Sk{12}.sk1...
    ;DS_Sk{13}.sk1;DS_Sk{14}.sk1;DS_Sk{15}.sk1;DS_Sk{16}.sk1];
SS = cat(2,SS,DSi_Sk);
SS = round(SS,3);

MaxCor = [DS_maxcor{1}.CCall(101);DS_maxcor{2}.CCall(101);DS_maxcor{3}.CCall(101);...
    DS_maxcor{4}.CCall(101);DS_maxcor{5}.CCall(101);DS_maxcor{6}.CCall(101);...
    DS_maxcor{7}.CCall(101);DS_maxcor{8}.CCall(101);DS_maxcor{9}.CCall(101);...
    DS_maxcor{10}.CCall(101);DS_maxcor{11}.CCall(101);DS_maxcor{12}.CCall(101)...
    ;DS_maxcor{13}.CCall(101);DS_maxcor{14}.CCall(101);DS_maxcor{15}.CCall(101)...
    ;DS_maxcor{16}.CCall(101)];
MaxCor = cat(2,MaxCor,DSi_maxcor);
MaxCor = round(MaxCor,3);

STD = [DS_STD{1};DS_STD{2};DS_STD{3};DS_STD{4};DS_STD{5};DS_STD{6};...
    DS_STD{7};DS_STD{8};DS_STD{9};DS_STD{10};DS_STD{11};DS_STD{12}...
    ;DS_STD{13};DS_STD{14};DS_STD{15};DS_STD{16}];
STD = cat(2,STD,DSi_std_M);
STD = round(STD,3);

DSpreadStats = table(Model,f_Band,RMSE,SS,MaxCor,STD)

clear RMSE SS MaxCor STD


% OBSERVED STD FOR INDIVIDUAL BUOYS

EnergySTD = Ei_std_O;
DirSTD = Di_std_O;
DSpreadSTD = DSi_std_O;

iObsSTD = table(Model,f_Band,EnergySTD,DirSTD,DSpreadSTD)

clear EnergySTD DirSTD DSpreadSTD


%% Save and Clear Vars


save('InterpModelObs.mat','-append','Buoy')


save('obsVSmod_STATS.mat','-append','EnergyStats','DirectionStats',...
    'DSpreadStats','E_OSTD','D_OSTD','DS_OSTD','iObsSTD')

saveas(figure(1),'Figures/ObsvsMod_intE.jpeg')
saveas(figure(2),'Figures/ObsvsMod_EWMdir.jpeg')
saveas(figure(3),'Figures/ObsvsMod_EWMspread.jpeg')



clear i MintE1 MintE2 MintE3 MintE4 OintE1 OintE2 OintE3 OintE4 dim str 
clear M_D_EWM1 M_D_EWM2 M_D_EWM3 M_D_EWM4 O_D_EWM1 O_D_EWM2 O_D_EWM3 O_D_EWM4
clear M_DS_EWM1 M_DS_EWM2 M_DS_EWM3 M_DS_EWM4 O_DS_EWM1 O_DS_EWM2 O_DS_EWM3 O_DS_EWM4
clear xL yL a Bnames SWAN_CATO SWAN_CATM SWANsm_CATO SWANsm_CATM COUP_CATM COUP_CATO
clear D_maxcor D_rmse D_Sk DS_maxcor DS_rmse DS_Sk E_maxcor E_rmse E_Sk
clear pformat one_to_one Model f_Band lightgrey transformat

clear D_OSTD D_STD Di_maxcor Di_rmse Di_Sk Di_std_M Di_std_O DS_OSTD
clear DS_STD DSi_maxcor DSi_rmse DSi_Sk DSi_std_M DSi_std_O E_OSTD E_STD
clear Ei_maxcor Ei_rmse Ei_Sk Ei_Sk4 Ei_std_M Ei_std_O fdivfp maxcor Mod
clear MSize Obs rmse Sk












