%% calc_MvsO_Stats_All.m
%
% Noah Clark            7/25/24
%
%
% Purpose: Calculate the means per hour over most of the spectra 
%          (containing the majority of the energy) for each of the models
%          and the observed. Then compare them by using the mod_error
%          function to calculate the statistics 
%
%
%

%% Notes

% Freq Range: 0.04 Hz - 0.23 Hz (ind: 1-20)


%% Preliminaries

clc;clear;

addpath('../')  % path to load data
addpath('../../Stat Functions') % the functions for the statistical analysis
addpath('C:\Users\nsc4750\OneDrive - UNC-Wilmington\CMS Work\Summer (0516-0728)\Main Functions') %path for TED function
load('InterpModelObs.mat')


Bnames = fields(Buoy.Obs); % to list buoy names

iF = 1:20; % indecies for frequency range (0.04 Hz - 0.23 Hz)




%% ----------------------------------------------------------------
%%              DETERMINE STATISTICS AT EACH BUOY 
%% ----------------------------------------------------------------


%% Calculate Averages for E, Dir, and Spread
% Energy - Integrated energy over frequency range
% Direction - EWM direction over frequency range
% Spread - EWM spread over frequency range


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
            
            % Integrated Energy
            intE(j) = sum(STODO.See(iF,j).*0.01);
            
            % EWM Direction
            m0 = intE(j);
            Dm1 = nansum(STODO.See(iF,j).*STODO.Dir(iF,j)).*0.01;
            D_EWM(j) = Dm1/m0;
            
            %EWM Spread
            DSm1 = nansum(STODO.See(iF,j).*STODO.Dspread(iF,j)).*0.01;
            DS_EWM(j) = DSm1/m0;

        end
        
        if MO == 1
            eval(sprintf('Buoy.Obs.%s.intE_All = intE;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMDir_All = D_EWM;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMDspread_All = DS_EWM;',Bnames{i}))
        elseif MO == 2
            eval(sprintf('Buoy.SWAN050.%s.intE_All = intE;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMDir_All = D_EWM;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMDspread_All = DS_EWM;',Bnames{i}))
        elseif MO == 3
            eval(sprintf('Buoy.SWAN050sm.%s.intE_All = intE;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMDir_All = D_EWM;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMDspread_All = DS_EWM;',Bnames{i}))
        elseif MO == 4
            eval(sprintf('Buoy.COUP050.%s.intE_All = intE;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMDir_All = D_EWM;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMDspread_All = DS_EWM;',Bnames{i}))
        else
            eval(sprintf('Buoy.SWAN075.%s.intE_All = intE;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMDir_All = D_EWM;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMDspread_All = DS_EWM;',Bnames{i}))
        end

    end
end

clear intE m0 Dm1 D_EWM DSm1 DS_EWM


%% Calculate Statistics for E, Dir, and Spread

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
        
        % Energy
        [~,Sk,maxcor,rmse] = mod_error(Mod.intE_All,Obs.intE_All);
        Ei_Sk(MO,j) = Sk.sk1;
        Ei_maxcor(MO,j) = maxcor.CCall(101);
        Ei_rmse(MO,j) = rmse.r1;
        Eistd_M(MO,j) = nanstd(Mod.intE_All);
        Eistd_O(MO,j) = nanstd(Obs.intE_All);
        
        % Direction
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMDir_All,Obs.EWMDir_All);
        Di_Sk(MO,j) = Sk.sk1;
        Di_maxcor(MO,j) = maxcor.CCall(101);
        Di_rmse(MO,j) = rmse.r1;
        Distd_M(MO,j) = nanstd(Mod.EWMDir_All);
        Distd_O(MO,j) = nanstd(Obs.EWMDir_All);
        
        % Spread
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMDspread_All,Obs.EWMDspread_All);
        DSi_Sk(MO,j) = Sk.sk1;
        DSi_maxcor(MO,j) = maxcor.CCall(101);
        DSi_rmse(MO,j) = rmse.r1;
        DSistd_M(MO,j) = nanstd(Mod.EWMDspread_All);
        DSistd_O(MO,j) = nanstd(Obs.EWMDspread_All);
        
    end
end







%% ----------------------------------------------------------------
%%              DETERMINE STATISTICS BETWEEN EACH BUOY 
%% ----------------------------------------------------------------


%% Calculate Averages for TED, Diff-Dir, and Diff-Spread
% All EWM

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
            % TED
            [TED,~,~,~] = NC_ObsDissDF(STODO1.See(iF,j),STODO2.See(iF,j),...
                STODO1.Dir(iF,j),STODO2.Dir(iF,j),Buoy.freq(iF),...
                STODO1.utm,STODO2.utm,STODO1.depth(j),STODO2.depth(j),0.01);
            avgSee = (STODO1.See(iF,j) + STODO2.See(iF,j))./2;
            Am0 = sum(avgSee).*0.01;
            m1TED = sum(avgSee.*TED).*0.01;
            EWM_TED(j) = m1TED./Am0;
            
            % Difference in Direction
            diffD = STODO2.Dir(iF,j) - STODO1.Dir(iF,j); 
            Dm1 = nansum(diffD.*avgSee).*0.01;
            D_EWM(j) = Dm1./Am0;
            
            % Difference in Spread
            diffDS = STODO2.Dspread(iF,j) - STODO1.Dspread(iF,j); 
            DSm1 = nansum(diffDS.*avgSee).*0.01;
            DS_EWM(j) = DSm1./Am0;
            
        end

        if MO == 1
            eval(sprintf('Buoy.Obs.%s.EWM_TED_All = EWM_TED;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMdiffDir_All = D_EWM;',Bnames{i}))
            eval(sprintf('Buoy.Obs.%s.EWMdiffDspread_All = DS_EWM;',Bnames{i}))
        elseif MO == 2
            eval(sprintf('Buoy.SWAN050.%s.EWM_TED_All = EWM_TED;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMdiffDir_All = D_EWM;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050.%s.EWMdiffDspread_All = DS_EWM;',Bnames{i}))
        elseif MO == 3
            eval(sprintf('Buoy.SWAN050sm.%s.EWM_TED_All = EWM_TED;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMdiffDir_All = D_EWM;',Bnames{i}))
            eval(sprintf('Buoy.SWAN050sm.%s.EWMdiffDspread_All = DS_EWM;',Bnames{i}))
        elseif MO == 4
            eval(sprintf('Buoy.COUP050.%s.EWM_TED_All = EWM_TED;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMdiffDir_All = D_EWM;',Bnames{i}))
            eval(sprintf('Buoy.COUP050.%s.EWMdiffDspread_All = DS_EWM;',Bnames{i}))
        else
            eval(sprintf('Buoy.SWAN075.%s.EWM_TED_All = EWM_TED;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMdiffDir_All = D_EWM;',Bnames{i}))
            eval(sprintf('Buoy.SWAN075.%s.EWMdiffDspread_All = DS_EWM;',Bnames{i}))
        end
        
    end
end


clear TED avgSee Am0 m1TED EWM_TED diffD Dm1 D_EWM diffDS DSm1 DS_EWM


%% Calculate Statistics for TED, Diff-Dir, and Diff-Spread


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
        
        % TED
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWM_TED_All,Obs.EWM_TED_All);
        TEDi_Sk(MO,j) = Sk.sk1;
        TEDi_maxcor(MO,j) = maxcor.CCall(101);
        TEDi_rmse(MO,j) = rmse.r1;
        TEDistd_M(MO,j) = nanstd(Mod.EWM_TED_All);
        TEDistd_O(MO,j) = nanstd(Obs.EWM_TED_All);
        
        % Difference in Direction
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMdiffDir_All,Obs.EWMdiffDir_All);
        dDi_Sk(MO,j) = Sk.sk1;
        dDi_maxcor(MO,j) = maxcor.CCall(101);
        dDi_rmse(MO,j) = rmse.r1;
        dDistd_M(MO,j) = nanstd(Mod.EWMdiffDir_All);
        dDistd_O(MO,j) = nanstd(Obs.EWMdiffDir_All);
        
        % Difference in Spread
        [~,Sk,maxcor,rmse] = mod_error(Mod.EWMdiffDspread_All,Obs.EWMdiffDspread_All);
        dDSi_Sk(MO,j) = Sk.sk1;
        dDSi_maxcor(MO,j) = maxcor.CCall(101);
        dDSi_rmse(MO,j) = rmse.r1;
        dDSistd_M(MO,j) = nanstd(Mod.EWMdiffDspread_All);
        dDSistd_O(MO,j) = nanstd(Obs.EWMdiffDspread_All);
        
    end
end




%% ----------------------------------------------------------------
%%          COMBINE ALL STATISTICS INTO A TABLE AND SAVE
%% ----------------------------------------------------------------


%% Create Tables

Model = ["SWAN050";"SWAN050sm";"COUP050";"SWAN075"];


    % ENERGY TABLE:
RMSE = Ei_rmse;
SS = Ei_Sk;
MaxCor = Ei_maxcor;
STD = Eistd_M;

Stats_FS.EnergyStats = table(Model,RMSE,SS,MaxCor,STD);

clear RMSE SS MaxCor STD


    %  DIRECTION TABLE:
RMSE = Di_rmse;
SS = Di_Sk;
MaxCor = Di_maxcor;
STD = Distd_M;

Stats_FS.DirectionStats = table(Model,RMSE,SS,MaxCor,STD);

clear RMSE SS MaxCor STD


    % SPREAD TABLE:
RMSE = DSi_rmse;
SS = DSi_Sk;
MaxCor = DSi_maxcor;
STD = DSistd_M;

Stats_FS.SpreadStats = table(Model,RMSE,SS,MaxCor,STD);

clear RMSE SS MaxCor STD


    % TED TABLE:
RMSE = TEDi_rmse;
SS = TEDi_Sk;
MaxCor = TEDi_maxcor;
STD = TEDistd_M;

Stats_FS.TEDStats = table(Model,RMSE,SS,MaxCor,STD);

clear RMSE SS MaxCor STD


    % DIFFERENCE IN DIRECTION TABLE:
RMSE = dDi_rmse;
SS = dDi_Sk;
MaxCor = dDi_maxcor;
STD = dDistd_M;

Stats_FS.diffDirStats = table(Model,RMSE,SS,MaxCor,STD);

clear RMSE SS MaxCor STD


    % DIFFERENCE IN SPREAD TABLE:
RMSE = dDSi_rmse;
SS = dDSi_Sk;
MaxCor = dDSi_maxcor;
STD = dDSistd_M;

Stats_FS.diffSpreadStats = table(Model,RMSE,SS,MaxCor,STD);

clear RMSE SS MaxCor STD


% TABLE FOR THE OBSERVED STD
E_STD = Eistd_O;
Dir_STD = Distd_O;
Spread_STD = DSistd_O;
TED_STD = TEDistd_O;
diffDir_STD = dDistd_O;
diffSpread_STD = dDSistd_O;

Stats_FS.iObsSTD = table(Model,E_STD,Dir_STD,Spread_STD,TED_STD,diffDir_STD,diffSpread_STD);



Stats_FS


%% Save and Clear Variables



save('Stats_FS.mat','Stats_FS')
























