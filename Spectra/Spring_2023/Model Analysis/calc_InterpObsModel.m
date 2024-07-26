%% calc_InterpObsModel.m
%
% Noah Clark         5/3/2024
%
%
% Purpose: Interpolate the observed and model data to the new frequency
%          vector. Also, cut the times to make up the two weeks ran for
%          SWAN050sm.
%
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
%
% Note:

% CR = 285
% X = 293



%% Preliminaries

clc;clear;

% Model Data
addpath('Model Results')
load('ModelResults.mat')

% Observed Buoy Data
addpath('../')
load('WBvariables.mat','BSee','XSee','Btime','Xtime','BEMEM','XEMEM',...
    'BDspread','Xdepth','Bdepth','Xutm','Butm')
% Observed ADCP Data
load('SM&ADCP_All.mat','ADCP')

% Get Spreads for Asilomar Buoys
addpath('../../Summer (0516-0728)/Start Data')
load('roxsi_spotter_L2_X01_1151_reduced.mat')
XDspread{1} = spotterL2.EMEM.meandirspread_f;
load('roxsi_spotter_L2_X03_1157_reduced.mat')
XDspread{2} = spotterL2.EMEM.meandirspread_f;
load('roxsi_spotter_L2_X04_1155_reduced.mat')
XDspread{3} = spotterL2.EMEM.meandirspread_f;

% INTERPOLATE EVERYTHING ONTO THIS FREQ VECTOR INSTEAD
ffNEW = 0.04:0.01:0.5;
ff = (1:129).*0.0098;

%% Import New Observed Directions and Spreads

addpath('ObsDirs&Spreads')
load('waves_allmoorings.mat')
NDDS_time = B01.time;
NDDS_time = datetime(NDDS_time,'ConvertFrom','datenum');
NDDS_f = B01.frequency(2:end);

BNormD = 285;
XNormD = 293;


%% Add NaNs in B01 blank space

B01See_1 = BSee{1}(:,1:397);
B01See_2 = BSee{1}(:,398:end);

B01Dir_1 = BEMEM{1}(:,1:397);
B01Dir_2 = BEMEM{1}(:,398:end);

B01Spread_1 = BDspread{1}(:,1:397);
B01Spread_2 = BDspread{1}(:,398:end);

B01NaNs = nan(129,107);


BSee{1} = cat(2,B01See_1,B01NaNs,B01See_2);
BEMEM{1} = cat(2,B01Dir_1,B01NaNs,B01Dir_2);
BDspread{1} = cat(2,B01Spread_1,B01NaNs,B01Spread_2);

clear B01See_1 B01See_2 B01Dir_1 B01Dir_2 B01Spread_1 B01Spread_2 B01NaNs


%% Cut times to all be the same
% Buoy Time: 15 June @12:00 - 20 July @03:00
% B01 Buy Time: 15 June @12:00 - 02 July @00:00 && 06 July @12:00 - 20 July @03:00
%       - changed to have nans in the blank space
% SWAN050 Time: 17 June @00:00 - 20 July @00:00
% SWAN050sm Time: 17 June @00:00 - 17 July @00:00
% COUP050 Time: 17 June @00:00 - 17 July @00:00
% ADCP Time: 21 June @19:00 - 25 July @10:00
% NDDS Time: 15 June @20:00 - 21 July @05:00

% CUT EVERYTHING: 21 June @19:00 - 17 July @00:00
Tvec = 1:606;
Oadd = 151;
ADCPadd = 0;
Madd = 115;
NDDSadd = 143;


    % Observed
Buoy.Obs.B01.See = BSee{1}(2:end,Tvec+Oadd);
Buoy.Obs.B03.See = BSee{2}(2:end,Tvec+Oadd);
Buoy.Obs.B05.See = BSee{3}(2:end,Tvec+Oadd);
Buoy.Obs.B10.See = ADCP.B10.See(:,Tvec+ADCPadd);
Buoy.Obs.B13.See = ADCP.B13.See(:,Tvec+ADCPadd);

% Buoy.Obs.B01.Dir = BEMEM{1}(2:end,Tvec+Oadd);
% Buoy.Obs.B03.Dir = BEMEM{2}(2:end,Tvec+Oadd);
% Buoy.Obs.B05.Dir = BEMEM{3}(2:end,Tvec+Oadd);
% Buoy.Obs.B10.Dir = ones(length(ffNEW),length(Tvec)) + NaN; % no ADCP directions
% Buoy.Obs.B13.Dir = ones(length(ffNEW),length(Tvec)) + NaN; % no ADCP directions
Buoy.Obs.B01.Dir = -B01.meandir2(Tvec+NDDSadd,2:end)' + BNormD;
Buoy.Obs.B03.Dir = -B03.meandir2(Tvec+NDDSadd,2:end)' + BNormD;
Buoy.Obs.B05.Dir = -B05.meandir2(Tvec+NDDSadd,2:end)' + BNormD;
Buoy.Obs.B10.Dir = -B10.meandir2(Tvec+NDDSadd,2:end)' + BNormD;
Buoy.Obs.B13.Dir = -B13.meandir2(Tvec+NDDSadd,2:end)' + BNormD;

% Buoy.Obs.B01.Dspread = BDspread{1}(2:end,Tvec+Oadd);
% Buoy.Obs.B03.Dspread = BDspread{2}(2:end,Tvec+Oadd);
% Buoy.Obs.B05.Dspread = BDspread{3}(2:end,Tvec+Oadd);
% Buoy.Obs.B10.Dspread = ones(length(ffNEW),length(Tvec)) + NaN; % no ADCP spreads
% Buoy.Obs.B13.Dspread = ones(length(ffNEW),length(Tvec)) + NaN; % no ADCP spreads
Buoy.Obs.B01.Dspread = B01.dirspread2(Tvec+NDDSadd,2:end)';
Buoy.Obs.B03.Dspread = B03.dirspread2(Tvec+NDDSadd,2:end)';
Buoy.Obs.B05.Dspread = B05.dirspread2(Tvec+NDDSadd,2:end)';
Buoy.Obs.B10.Dspread = B10.dirspread2(Tvec+NDDSadd,2:end)';
Buoy.Obs.B13.Dspread = B13.dirspread2(Tvec+NDDSadd,2:end)';

Buoy.Obs.X01.See = XSee{1}(2:end,Tvec+Oadd);
Buoy.Obs.X03.See = XSee{2}(2:end,Tvec+Oadd);
Buoy.Obs.X04.See = XSee{3}(2:end,Tvec+Oadd);
Buoy.Obs.X05.See = ADCP.X05.See(:,Tvec+ADCPadd);

% Buoy.Obs.X01.Dir = XEMEM{1}(2:end,Tvec+Oadd);
% Buoy.Obs.X03.Dir = XEMEM{2}(2:end,Tvec+Oadd);
% Buoy.Obs.X04.Dir = XEMEM{3}(2:end,Tvec+Oadd);
% Buoy.Obs.X05.Dir = ones(length(ffNEW),length(Tvec)) + NaN; % no ADCP directions
Buoy.Obs.X01.Dir = -X01.meandir2(Tvec+NDDSadd,2:end)' + XNormD;
Buoy.Obs.X03.Dir = -X03.meandir2(Tvec+NDDSadd,2:end)' + XNormD;
Buoy.Obs.X04.Dir = -X04.meandir2(Tvec+NDDSadd,2:end)' + XNormD;
Buoy.Obs.X05.Dir = -X05.meandir2(Tvec+NDDSadd,2:end)' + XNormD;

% Buoy.Obs.X01.Dspread = XDspread{1}(2:end,Tvec+Oadd);
% Buoy.Obs.X03.Dspread = XDspread{2}(2:end,Tvec+Oadd);
% Buoy.Obs.X04.Dspread = XDspread{3}(2:end,Tvec+Oadd);
% Buoy.Obs.X05.Dspread = ones(length(ffNEW),length(Tvec)) + NaN; % no ADCP spreads
Buoy.Obs.X01.Dspread = X01.dirspread2(Tvec+NDDSadd,2:end)';
Buoy.Obs.X03.Dspread = X03.dirspread2(Tvec+NDDSadd,2:end)';
Buoy.Obs.X04.Dspread = X04.dirspread2(Tvec+NDDSadd,2:end)';
Buoy.Obs.X05.Dspread = X05.dirspread2(Tvec+NDDSadd,2:end)';

    % SWAN050
Buoy.SWAN050.B01.See = SWAN050.B01.WS(1:40,Tvec+Madd);
Buoy.SWAN050.B03.See = SWAN050.B03.WS(1:40,Tvec+Madd);
Buoy.SWAN050.B05.See = SWAN050.B05.WS(1:40,Tvec+Madd);
Buoy.SWAN050.B10.See = SWAN050.B10.WS(1:40,Tvec+Madd);
Buoy.SWAN050.B13.See = SWAN050.B13.WS(1:40,Tvec+Madd);

Buoy.SWAN050.B01.Dir = SWAN050.B01.TTm(1:40,Tvec+Madd);
Buoy.SWAN050.B03.Dir = SWAN050.B03.TTm(1:40,Tvec+Madd);
Buoy.SWAN050.B05.Dir = SWAN050.B05.TTm(1:40,Tvec+Madd);
Buoy.SWAN050.B10.Dir = SWAN050.B10.TTm(1:40,Tvec+Madd);
Buoy.SWAN050.B13.Dir = SWAN050.B13.TTm(1:40,Tvec+Madd);

Buoy.SWAN050.B01.Dspread = SWAN050.B01.Spread(1:40,Tvec+Madd);
Buoy.SWAN050.B03.Dspread = SWAN050.B03.Spread(1:40,Tvec+Madd);
Buoy.SWAN050.B05.Dspread = SWAN050.B05.Spread(1:40,Tvec+Madd);
Buoy.SWAN050.B10.Dspread = SWAN050.B10.Spread(1:40,Tvec+Madd);
Buoy.SWAN050.B13.Dspread = SWAN050.B13.Spread(1:40,Tvec+Madd);

Buoy.SWAN050.X01.See = SWAN050.X01.WS(1:40,Tvec+Madd);
Buoy.SWAN050.X03.See = SWAN050.X03.WS(1:40,Tvec+Madd);
Buoy.SWAN050.X04.See = SWAN050.X04.WS(1:40,Tvec+Madd);
Buoy.SWAN050.X05.See = SWAN050.X05.WS(1:40,Tvec+Madd);

Buoy.SWAN050.X01.Dir = SWAN050.X01.TTm(1:40,Tvec+Madd);
Buoy.SWAN050.X03.Dir = SWAN050.X03.TTm(1:40,Tvec+Madd);
Buoy.SWAN050.X04.Dir = SWAN050.X04.TTm(1:40,Tvec+Madd);
Buoy.SWAN050.X05.Dir = SWAN050.X05.TTm(1:40,Tvec+Madd);

Buoy.SWAN050.X01.Dspread = SWAN050.X01.Spread(1:40,Tvec+Madd);
Buoy.SWAN050.X03.Dspread = SWAN050.X03.Spread(1:40,Tvec+Madd);
Buoy.SWAN050.X04.Dspread = SWAN050.X04.Spread(1:40,Tvec+Madd);
Buoy.SWAN050.X05.Dspread = SWAN050.X05.Spread(1:40,Tvec+Madd);

    % SWAN050sm
Buoy.SWAN050sm.B01.See = SWAN050sm.B01.WS(1:40,Tvec+Madd);
Buoy.SWAN050sm.B03.See = SWAN050sm.B03.WS(1:40,Tvec+Madd);
Buoy.SWAN050sm.B05.See = SWAN050sm.B05.WS(1:40,Tvec+Madd);
Buoy.SWAN050sm.B10.See = SWAN050sm.B10.WS(1:40,Tvec+Madd);
Buoy.SWAN050sm.B13.See = SWAN050sm.B13.WS(1:40,Tvec+Madd);

Buoy.SWAN050sm.B01.Dir = SWAN050sm.B01.TTm(1:40,Tvec+Madd);
Buoy.SWAN050sm.B03.Dir = SWAN050sm.B03.TTm(1:40,Tvec+Madd);
Buoy.SWAN050sm.B05.Dir = SWAN050sm.B05.TTm(1:40,Tvec+Madd);
Buoy.SWAN050sm.B10.Dir = SWAN050sm.B10.TTm(1:40,Tvec+Madd);
Buoy.SWAN050sm.B13.Dir = SWAN050sm.B13.TTm(1:40,Tvec+Madd);

Buoy.SWAN050sm.B01.Dspread = SWAN050sm.B01.Spread(1:40,Tvec+Madd);
Buoy.SWAN050sm.B03.Dspread = SWAN050sm.B03.Spread(1:40,Tvec+Madd);
Buoy.SWAN050sm.B05.Dspread = SWAN050sm.B05.Spread(1:40,Tvec+Madd);
Buoy.SWAN050sm.B10.Dspread = SWAN050sm.B10.Spread(1:40,Tvec+Madd);
Buoy.SWAN050sm.B13.Dspread = SWAN050sm.B13.Spread(1:40,Tvec+Madd);

Buoy.SWAN050sm.X01.See = SWAN050sm.X01.WS(1:40,Tvec+Madd);
Buoy.SWAN050sm.X03.See = SWAN050sm.X03.WS(1:40,Tvec+Madd);
Buoy.SWAN050sm.X04.See = SWAN050sm.X04.WS(1:40,Tvec+Madd);
Buoy.SWAN050sm.X05.See = SWAN050sm.X05.WS(1:40,Tvec+Madd);

Buoy.SWAN050sm.X01.Dir = SWAN050sm.X01.TTm(1:40,Tvec+Madd);
Buoy.SWAN050sm.X03.Dir = SWAN050sm.X03.TTm(1:40,Tvec+Madd);
Buoy.SWAN050sm.X04.Dir = SWAN050sm.X04.TTm(1:40,Tvec+Madd);
Buoy.SWAN050sm.X05.Dir = SWAN050sm.X05.TTm(1:40,Tvec+Madd);

Buoy.SWAN050sm.X01.Dspread = SWAN050sm.X01.Spread(1:40,Tvec+Madd);
Buoy.SWAN050sm.X03.Dspread = SWAN050sm.X03.Spread(1:40,Tvec+Madd);
Buoy.SWAN050sm.X04.Dspread = SWAN050sm.X04.Spread(1:40,Tvec+Madd);
Buoy.SWAN050sm.X05.Dspread = SWAN050sm.X05.Spread(1:40,Tvec+Madd);

    % COUP050
Buoy.COUP050.B01.See = COUP050.B01.WS(1:40,Tvec+Madd);
Buoy.COUP050.B03.See = COUP050.B03.WS(1:40,Tvec+Madd);
Buoy.COUP050.B05.See = COUP050.B05.WS(1:40,Tvec+Madd);
Buoy.COUP050.B10.See = COUP050.B10.WS(1:40,Tvec+Madd);
Buoy.COUP050.B13.See = COUP050.B13.WS(1:40,Tvec+Madd);

Buoy.COUP050.B01.Dir = COUP050.B01.TTm(1:40,Tvec+Madd);
Buoy.COUP050.B03.Dir = COUP050.B03.TTm(1:40,Tvec+Madd);
Buoy.COUP050.B05.Dir = COUP050.B05.TTm(1:40,Tvec+Madd);
Buoy.COUP050.B10.Dir = COUP050.B10.TTm(1:40,Tvec+Madd);
Buoy.COUP050.B13.Dir = COUP050.B13.TTm(1:40,Tvec+Madd);

Buoy.COUP050.B01.Dspread = COUP050.B01.Spread(1:40,Tvec+Madd);
Buoy.COUP050.B03.Dspread = COUP050.B03.Spread(1:40,Tvec+Madd);
Buoy.COUP050.B05.Dspread = COUP050.B05.Spread(1:40,Tvec+Madd);
Buoy.COUP050.B10.Dspread = COUP050.B10.Spread(1:40,Tvec+Madd);
Buoy.COUP050.B13.Dspread = COUP050.B13.Spread(1:40,Tvec+Madd);

Buoy.COUP050.X01.See = COUP050.X01.WS(1:40,Tvec+Madd);
Buoy.COUP050.X03.See = COUP050.X03.WS(1:40,Tvec+Madd);
Buoy.COUP050.X04.See = COUP050.X04.WS(1:40,Tvec+Madd);
Buoy.COUP050.X05.See = COUP050.X05.WS(1:40,Tvec+Madd);

Buoy.COUP050.X01.Dir = COUP050.X01.TTm(1:40,Tvec+Madd);
Buoy.COUP050.X03.Dir = COUP050.X03.TTm(1:40,Tvec+Madd);
Buoy.COUP050.X04.Dir = COUP050.X04.TTm(1:40,Tvec+Madd);
Buoy.COUP050.X05.Dir = COUP050.X05.TTm(1:40,Tvec+Madd);

Buoy.COUP050.X01.Dspread = COUP050.X01.Spread(1:40,Tvec+Madd);
Buoy.COUP050.X03.Dspread = COUP050.X03.Spread(1:40,Tvec+Madd);
Buoy.COUP050.X04.Dspread = COUP050.X04.Spread(1:40,Tvec+Madd);
Buoy.COUP050.X05.Dspread = COUP050.X05.Spread(1:40,Tvec+Madd);

    % SWAN075
Buoy.SWAN075.B01.See = SWAN075.B01.WS(1:40,Tvec+Madd);
Buoy.SWAN075.B03.See = SWAN075.B03.WS(1:40,Tvec+Madd);
Buoy.SWAN075.B05.See = SWAN075.B05.WS(1:40,Tvec+Madd);
Buoy.SWAN075.B10.See = SWAN075.B10.WS(1:40,Tvec+Madd);
Buoy.SWAN075.B13.See = SWAN075.B13.WS(1:40,Tvec+Madd);

Buoy.SWAN075.B01.Dir = SWAN075.B01.TTm(1:40,Tvec+Madd);
Buoy.SWAN075.B03.Dir = SWAN075.B03.TTm(1:40,Tvec+Madd);
Buoy.SWAN075.B05.Dir = SWAN075.B05.TTm(1:40,Tvec+Madd);
Buoy.SWAN075.B10.Dir = SWAN075.B10.TTm(1:40,Tvec+Madd);
Buoy.SWAN075.B13.Dir = SWAN075.B13.TTm(1:40,Tvec+Madd);

Buoy.SWAN075.B01.Dspread = SWAN075.B01.Spread(1:40,Tvec+Madd);
Buoy.SWAN075.B03.Dspread = SWAN075.B03.Spread(1:40,Tvec+Madd);
Buoy.SWAN075.B05.Dspread = SWAN075.B05.Spread(1:40,Tvec+Madd);
Buoy.SWAN075.B10.Dspread = SWAN075.B10.Spread(1:40,Tvec+Madd);
Buoy.SWAN075.B13.Dspread = SWAN075.B13.Spread(1:40,Tvec+Madd);

Buoy.SWAN075.X01.See = SWAN075.X01.WS(1:40,Tvec+Madd);
Buoy.SWAN075.X03.See = SWAN075.X03.WS(1:40,Tvec+Madd);
Buoy.SWAN075.X04.See = SWAN075.X04.WS(1:40,Tvec+Madd);
Buoy.SWAN075.X05.See = SWAN075.X05.WS(1:40,Tvec+Madd);

Buoy.SWAN075.X01.Dir = SWAN075.X01.TTm(1:40,Tvec+Madd);
Buoy.SWAN075.X03.Dir = SWAN075.X03.TTm(1:40,Tvec+Madd);
Buoy.SWAN075.X04.Dir = SWAN075.X04.TTm(1:40,Tvec+Madd);
Buoy.SWAN075.X05.Dir = SWAN075.X05.TTm(1:40,Tvec+Madd);

Buoy.SWAN075.X01.Dspread = SWAN075.X01.Spread(1:40,Tvec+Madd);
Buoy.SWAN075.X03.Dspread = SWAN075.X03.Spread(1:40,Tvec+Madd);
Buoy.SWAN075.X04.Dspread = SWAN075.X04.Spread(1:40,Tvec+Madd);
Buoy.SWAN075.X05.Dspread = SWAN075.X05.Spread(1:40,Tvec+Madd);


% Master Time Vector
Buoy.time = Btime{2}(Tvec+Oadd); % 21 June @19:00 - 01 July @00:00


%% Interpolate
% - Interpolate the observed and model data to the new frequency vector
%    (freq: 0.04 - 0.50 Hz)


Obsff = 0.0098:0.0098:(0.0098*128);
Obsff = Obsff';
ADCPff = ADCP.freq';
Modelff = SWAN050.B01.f(1:40)';

Buoy.freq = ffNEW';


Bnames = fields(SWAN050);

    % Observed
for i = 1:9
    
    if i~=4 && i~=5 && i~=9 %assign the correct original freq vector
        Cff = Obsff;
    else
        Cff = ADCPff;
    end
    
    eval(sprintf('STODO = Buoy.Obs.%s;',Bnames{i}))
    
    intSee = interp1(Cff,STODO.See,Buoy.freq);
    eval(sprintf('Buoy.Obs.%s.See = intSee;',Bnames{i}))
    
    intDir = interp1(NDDS_f,STODO.Dir,Buoy.freq);
    eval(sprintf('Buoy.Obs.%s.Dir = intDir;',Bnames{i}))
    intSpread = interp1(NDDS_f,STODO.Dspread,Buoy.freq);
    eval(sprintf('Buoy.Obs.%s.Dspread = intSpread;',Bnames{i}))

%     if i~=4 && i~=5 && i~=9
%         intDir = interp1(Cff,STODO.Dir,Buoy.freq);
%         eval(sprintf('Buoy.Obs.%s.Dir = intDir;',Bnames{i}))
%         intSpread = interp1(Cff,STODO.Dspread,Buoy.freq);
%         eval(sprintf('Buoy.Obs.%s.Dspread = intSpread;',Bnames{i}))
%     end
    
end


    % Modeled
for i = 1:9
    
        % SWAN050
    eval(sprintf('STODO = Buoy.SWAN050.%s;',Bnames{i}))
    intSee = interp1(Modelff,STODO.See,Buoy.freq);
    eval(sprintf('Buoy.SWAN050.%s.See = intSee;',Bnames{i}))
    intDir = interp1(Modelff,STODO.Dir,Buoy.freq);
    eval(sprintf('Buoy.SWAN050.%s.Dir = intDir;',Bnames{i}))
    intSpread = interp1(Modelff,STODO.Dspread,Buoy.freq);
    eval(sprintf('Buoy.SWAN050.%s.Dspread = intSpread;',Bnames{i}))
    
            % SWAN050sm
    eval(sprintf('STODO = Buoy.SWAN050sm.%s;',Bnames{i}))
    intSee = interp1(Modelff,STODO.See,Buoy.freq);
    eval(sprintf('Buoy.SWAN050sm.%s.See = intSee;',Bnames{i}))
    intDir = interp1(Modelff,STODO.Dir,Buoy.freq);
    eval(sprintf('Buoy.SWAN050sm.%s.Dir = intDir;',Bnames{i}))
    intSpread = interp1(Modelff,STODO.Dspread,Buoy.freq);
    eval(sprintf('Buoy.SWAN050sm.%s.Dspread = intSpread;',Bnames{i}))
    
            % COUP050
    eval(sprintf('STODO = Buoy.COUP050.%s;',Bnames{i}))
    intSee = interp1(Modelff,STODO.See,Buoy.freq);
    eval(sprintf('Buoy.COUP050.%s.See = intSee;',Bnames{i}))
    intDir = interp1(Modelff,STODO.Dir,Buoy.freq);
    eval(sprintf('Buoy.COUP050.%s.Dir = intDir;',Bnames{i}))
    intSpread = interp1(Modelff,STODO.Dspread,Buoy.freq);
    eval(sprintf('Buoy.COUP050.%s.Dspread = intSpread;',Bnames{i}))
    
            % SWAN075
    eval(sprintf('STODO = Buoy.SWAN075.%s;',Bnames{i}))
    intSee = interp1(Modelff,STODO.See,Buoy.freq);
    eval(sprintf('Buoy.SWAN075.%s.See = intSee;',Bnames{i}))
    intDir = interp1(Modelff,STODO.Dir,Buoy.freq);
    eval(sprintf('Buoy.SWAN075.%s.Dir = intDir;',Bnames{i}))
    intSpread = interp1(Modelff,STODO.Dspread,Buoy.freq);
    eval(sprintf('Buoy.SWAN075.%s.Dspread = intSpread;',Bnames{i}))
end



%% Add Depths and utm to the "Buoy" Variable

% Add Nans in missing section of B01
B01h1 = Bdepth{1}(1:397);
B01h2 = Bdepth{1}(398:end);
B01NaNs = nan(107,1);
Bdepth{1} = cat(1,B01h1,B01NaNs,B01h2);
clear B01h1 B01h2 B01NaNs

% Assign ADCP depth values to new structure to index
B_ADCP_h{1} = ADCP.B10.depth;
B_ADCP_h{2} = ADCP.B13.depth;
X_ADCP_h{1} = ADCP.X05.depth;
X_ADCP_h{2} = ADCP.X11.depth; % goes unused

% Assign ADCP utm values to new structure to index
B_ADCP_utm{1} = ADCP.B10.utm;
B_ADCP_utm{2} = ADCP.B13.utm;
X_ADCP_utm{1} = ADCP.X05.utm;
X_ADCP_utm{2} = ADCP.X11.utm; % goes unused


for i = 1:9
    if i == 1 || i == 2 || i == 3
        h = Bdepth{i}(Tvec+Oadd);
        utm = Butm{i};
    elseif i == 4 || i == 5
        h = B_ADCP_h{i-3}(Tvec+ADCPadd);
        utm = B_ADCP_utm{i-3};
    elseif i == 6 || i == 7 || i == 8
        h = Xdepth{i-5}(Tvec+Oadd);
        utm = Xutm{i-5};
    else
        h = X_ADCP_h{i-8}(Tvec+ADCPadd);
        utm = X_ADCP_utm{i-8};
    end
    
    eval(sprintf('Buoy.Obs.%s.depth = h;',Bnames{i}))
    eval(sprintf('Buoy.SWAN050.%s.depth = h;',Bnames{i}))
    eval(sprintf('Buoy.SWAN050sm.%s.depth = h;',Bnames{i}))
    eval(sprintf('Buoy.COUP050.%s.depth = h;',Bnames{i}))
    eval(sprintf('Buoy.SWAN075.%s.depth = h;',Bnames{i}))
    
    eval(sprintf('Buoy.Obs.%s.utm = utm;',Bnames{i}))
    eval(sprintf('Buoy.SWAN050.%s.utm = utm;',Bnames{i}))
    eval(sprintf('Buoy.SWAN050sm.%s.utm = utm;',Bnames{i}))
    eval(sprintf('Buoy.COUP050.%s.utm = utm;',Bnames{i}))
    eval(sprintf('Buoy.SWAN075.%s.utm = utm;',Bnames{i}))
end 


%% Clear and Save 

clear spotterL2 i utm h intDir intSee intSpread STODO Cff
clear B01 B03 B05 B10 B07 B13 X01 X03 X04 X05 X08

save('InterpModelObs.mat','Buoy')










