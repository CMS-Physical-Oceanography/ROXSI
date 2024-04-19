%% plot_DailyObsTEDs.m
%
% Noah Clark        4/~/2024
%
%
% Purpose: To perform observed TED analysis and determine the EWM per hour
%          in the sea and swell bands of these TED values. Then plot these
%          observed TED daily averages. This is done for all spotter buoys
%          and ADCPs at China Rock and Asilomar.
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


%% Preliminaries


clc;clear;

% Load data for China Rock
load('WBvariables.mat','BSee','Btime','Bdepth','BEMEM','Butm','BNormWaveDir')
% Load data for Asilomar
load('WBvariables.mat','XSee','Xtime','Xdepth','XEMEM','Xutm','XNormWaveDir')
% Load data for ADCPs
load('SM&ADCP_All.mat','ADCP')

% Path for Nielsen TED functions
addpath('C:\Users\nsc4750\OneDrive - UNC-Wilmington\CMS Work\Summer (0516-0728)\Main Functions')
%addpath('/Users/noahclark/Library/CloudStorage/OneDrive-UNC-Wilmington/CMS Work/Summer (0516-0728)/Main Functions')

    % Define Variables
ff = 0.0098:0.0098:(0.0098*128);
ffa = ff(1:35); %for adcp

Tday = datetime(2022,6,16):datetime(2022,7,19);
TdayB01 = [datetime(2022,6,16):datetime(2022,6,30),...
    datetime(2022,7,7):datetime(2022,7,19)];
TdayADCP = datetime(2022,6,25):datetime(2022,7,19);

ADCP.B10.depth(1:61) = ADCP.B13.depth(1:61) + 4.75; % to fill the NaN values


for i = 1:3
    BSee{i} = BSee{i}(2:end,:);
    XSee{i} = XSee{i}(2:end,:);
    
    BEMEM{i} = BEMEM{i}(2:end,:);
    XEMEM{i} = XEMEM{i}(2:end,:);
end


    % Interpolate ADCP Sees to spotter frequencies
ADCP.B10.See = interp1(ADCP.freq',ADCP.B10.See,ff');
ADCP.B10.See = ADCP.B10.See(1:35,66:666); %take out the times to fit with data from spotters (June 25 - July 19)
ADCP.B13.See = interp1(ADCP.freq',ADCP.B13.See,ff');
ADCP.B13.See = ADCP.B13.See(1:35,66:666);

ADCP.X05.See = interp1(ADCP.freq',ADCP.X05.See,ff');
ADCP.X05.See = ADCP.X05.See(1:35,66:666);
ADCP.X11.See = interp1(ADCP.freq',ADCP.X11.See,ff');
ADCP.X11.See = ADCP.X11.See(1:35,66:666);

% Assume the ADCP directions are shore normal
ADCP.B10.EMEM = ones(35,601).*BNormWaveDir;
ADCP.B13.EMEM = ADCP.B10.EMEM;
ADCP.X05.EMEM = ones(35,601).*XNormWaveDir;
ADCP.X11.EMEM = ADCP.X05.EMEM;

% Cut the ADCP depths to match the time of ADCP sees
ADCP.B10.depth = ADCP.B10.depth(66:666);
ADCP.B13.depth = ADCP.B13.depth(66:666);
ADCP.X05.depth = ADCP.X05.depth(66:666);
ADCP.X11.depth = ADCP.X11.depth(66:666);


%% China Rock

    % Determine all of the observed TEDs:

    % B01 to B03
[B.OTED{1},~,~,~] = NC_ObsDiss(BSee{1},BSee{2}(:,[1:397 505:end]),...
    BEMEM{1},BEMEM{2}(:,[1:397 505:end]),ff,Butm{1},Butm{2},...
    Bdepth{1},Bdepth{2}([1:397 505:end]));

    % B03 to B05
[B.OTED{2},~,~,~] = NC_ObsDiss(BSee{2},BSee{3},BEMEM{2},BEMEM{3},ff,...
    Butm{2},Butm{3},Bdepth{2},Bdepth{3});

    % B05 to B10
[B.OTED{3},~,~,~] = NC_ObsDiss(BSee{3}(1:35,217:817),ADCP.B10.See,...
    BEMEM{3}(1:35,217:817),ADCP.B10.EMEM,ff,Butm{3},ADCP.B10.utm,...
    Bdepth{3}(217:817),ADCP.B10.depth);

    % B10 to B13
[B.OTED{4},~,~,~] = NC_ObsDiss(ADCP.B10.See,ADCP.B13.See,ADCP.B10.EMEM,...
    ADCP.B13.EMEM,ff,ADCP.B10.utm,ADCP.B13.utm,ADCP.B10.depth,...
    ADCP.B13.depth);


        % AVERAGE THE SEES B/T BUOYS:
B.ASee{1} = (BSee{1} + BSee{2}(:,[1:397 505:end]))./2;
B.ASee{2} = (BSee{2} + BSee{3})./2;
B.ASee{3} = (BSee{3}(1:35,217:817) + ADCP.B10.See)./2;
B.ASee{4} = (ADCP.B10.See + ADCP.B13.See)./2;


        % EWM PER DAY FOR SEA AND SWELL:
for j = 1:4
    if j == 1 % B01 to B03
        N = 725;
    elseif j == 2 % B03 to B05
        N = 832;
    else % B05 to B10 & B10 to B13
        N = 601;
    end
    
    for i = 1:N
            % Sea
        m0Sea = trapz(ff(9:24),B.ASee{j}(9:24,i),1);
        m1Sea = trapz(ff(9:24),B.ASee{j}(9:24,i).*B.OTED{j}(9:24,i));
        B.EWM_OTED_Sea{j}(i) = m1Sea./m0Sea;
            % Swell
        m0Swell = trapz(ff(4:9),B.ASee{j}(4:9,i),1);
        m1Swell = trapz(ff(4:9),B.ASee{j}(4:9,i).*B.OTED{j}(4:9,i));
        B.EWM_OTED_Swell{j}(i) = m1Swell./m0Swell;
        
    end
end



        % AVERAGE PER DAY
% B01 to B03
% (**1st Section**)
for i = 1:15
    EWM_DD_OTED_Sea1(i) = mean(B.EWM_OTED_Sea{1}(i*24-23:i*24));
    EWM_DD_OTED_Swell1(i) = mean(B.EWM_OTED_Swell{1}(i*24-23:i*24));
end
% (**2nd Section**)
for i = 17.625:29.625
    EWM_DD_OTED_Sea2(i-16.625) = mean(B.EWM_OTED_Sea{1}(i*24-23:i*24));
    EWM_DD_OTED_Swell2(i-16.625) = mean(B.EWM_OTED_Swell{1}(i*24-23:i*24));
end

B.EWM_DD_OTEDSea{1} = cat(2,EWM_DD_OTED_Sea1,EWM_DD_OTED_Sea2);
B.EWM_DD_OTEDSwell{1} = cat(2,EWM_DD_OTED_Swell1,EWM_DD_OTED_Swell2);
clear EWM_DD_OTED_Sea1 EWM_DD_OTED_Sea2
clear EWM_DD_OTED_Swell1 EWM_DD_OTED_Swell2

% B03 to B05, B05 to B10, & B10 to B13
for j = 2:4
    if j == 2 % B03 to B05
        N = 34;
    else % B05 to B10 & B10 to B13
        N = 25;
    end
    
    for i = 1:N
        B.EWM_DD_OTEDSea{j}(i) = mean(B.EWM_OTED_Sea{j}(i*24-23:i*24));
        B.EWM_DD_OTEDSwell{j}(i) = mean(B.EWM_OTED_Swell{j}(i*24-23:i*24));
    end
end





%% Asilomar

% X01 to X03
[X.OTED{1},~,~,~] = NC_ObsDiss(XSee{1},XSee{2},XEMEM{1},XEMEM{2},ff,...
    Xutm{1},Xutm{2},Xdepth{1},Xdepth{2});

% X03 to X04
[X.OTED{2},~,~,~] = NC_ObsDiss(XSee{2},XSee{3},XEMEM{2},XEMEM{3},ff,...
    Xutm{2},Xutm{3},Xdepth{2},Xdepth{3});

% X04 to X05
[X.OTED{3},~,~,~] = NC_ObsDiss(XSee{3}(1:35,217:817),ADCP.X05.See,...
    XEMEM{3}(1:35,217:817),ADCP.X05.EMEM,ff,Xutm{3},ADCP.X05.utm,...
    Xdepth{3}(217:817),ADCP.X05.depth);

% X05 to X11
[X.OTED{4},~,~,~] = NC_ObsDiss(ADCP.X05.See,ADCP.X11.See,ADCP.X05.EMEM,...
    ADCP.X11.EMEM,ff,ADCP.X05.utm,ADCP.X11.utm,ADCP.X05.depth,...
    ADCP.X11.depth);


        % AVERAGE THE SEES B/T BUOYS:
X.ASee{1} = (XSee{1} + XSee{2})./2;
X.ASee{2} = (XSee{2} + XSee{3})./2;
X.ASee{3} = (XSee{3}(1:35,217:817) + ADCP.X05.See)./2;
X.ASee{4} = (ADCP.X05.See + ADCP.X11.See)./2;



        % EWM PER DAY FOR SEA AND SWELL:
for j = 1:4
    if j == 1 || j == 2 % X01 to X03 & X03 to X04
        N = 832;
    else % X04 to X05 & X05 to X11
        N = 601;
    end
    
    for i = 1:N
            % Sea
        m0Sea = trapz(ff(9:24),X.ASee{j}(9:24,i),1);
        m1Sea = trapz(ff(9:24),X.ASee{j}(9:24,i).*X.OTED{j}(9:24,i));
        X.EWM_OTED_Sea{j}(i) = m1Sea./m0Sea;
            % Swell
        m0Swell = trapz(ff(4:9),X.ASee{j}(4:9,i),1);
        m1Swell = trapz(ff(4:9),X.ASee{j}(4:9,i).*X.OTED{j}(4:9,i));
        X.EWM_OTED_Swell{j}(i) = m1Swell./m0Swell;
        
    end
end


        % AVERAGE PER DAY
for j = 1:4
    if j == 1 || j == 2 % X01 to X03 & X03 to X04
        N = 34;
    else % X04 to X05 & X05 to X11
        N = 25;
    end
    
    for i = 1:N
        X.EWM_DD_OTEDSea{j}(i) = nanmean(X.EWM_OTED_Sea{j}(i*24-23:i*24));
        X.EWM_DD_OTEDSwell{j}(i) = nanmean(X.EWM_OTED_Swell{j}(i*24-23:i*24));
    end
end


%% Average Observed TED b/t Buoys

for i = 1:4
    % China Rock
    B.Avg_OTED_Sea{i} = nanmean(B.EWM_DD_OTEDSea{i});
    B.Avg_OTED_Swell{i} = nanmean(B.EWM_DD_OTEDSwell{i});
    % Asilomar
    X.Avg_OTED_Sea{i} = nanmean(X.EWM_DD_OTEDSea{i});
    X.Avg_OTED_Swell{i} = nanmean(X.EWM_DD_OTEDSwell{i});
end



%% Plotting

figure(1);clf;

    % China Rock - Swell
subplot(2,2,1)
plot(TdayB01,B.EWM_DD_OTEDSwell{1},'bx','MarkerSize',10)
hold on; grid on;
plot(Tday,B.EWM_DD_OTEDSwell{2},'mx','MarkerSize',10)
plot(TdayADCP,B.EWM_DD_OTEDSwell{3},'gx','MarkerSize',10)
plot(TdayADCP,B.EWM_DD_OTEDSwell{4},'rx','MarkerSize',10)

title('China Rock: Swell Observed Energy Dissipation','fontsize',17)
legend('B01 - B03 (AVG=-0.216)','B03 - B05 (AVG=-0.011)',...
    'B05 - B10 (AVG=0.334)','B10 - B13 (AVG=2.483)','fontsize',10,...
    'location','north')
xlabel('Date','fontsize',15)
ylabel('Observed TED [W/m]','fontsize',15)
ylim([-3 20])
xlim([datetime(2022,6,15) datetime(2022,7,20)])

    % China Rock - Sea
subplot(2,2,2)
plot(TdayB01,B.EWM_DD_OTEDSea{1},'bx','MarkerSize',10)
hold on; grid on;
plot(Tday,B.EWM_DD_OTEDSea{2},'mx','MarkerSize',10)
plot(TdayADCP,B.EWM_DD_OTEDSea{3},'gx','MarkerSize',10)
plot(TdayADCP,B.EWM_DD_OTEDSea{4},'rx','MarkerSize',10)

title('China Rock: Sea Observed Energy Dissipation','fontsize',17)
legend('B01 - B03 (AVG=0.120)','B03 - B05 (AVG=0.487)',...
    'B05 - B10 (AVG=-0.393)','B10 - B13 (AVG=5.993)','fontsize',10,...
    'location','north')
xlabel('Date','fontsize',15)
ylabel('Observed TED [W/m]','fontsize',15)
ylim([-3 20])
xlim([datetime(2022,6,15) datetime(2022,7,20)])



    % Asilomar - Swell
subplot(2,2,3)
plot(Tday,X.EWM_DD_OTEDSwell{1},'bx','MarkerSize',10)
hold on; grid on;
plot(Tday,X.EWM_DD_OTEDSwell{2},'mx','MarkerSize',10)
plot(TdayADCP,X.EWM_DD_OTEDSwell{3},'gx','MarkerSize',10)
plot(TdayADCP,X.EWM_DD_OTEDSwell{4},'rx','MarkerSize',10)
title('Asilomar: Swell Observed Energy Dissipation','fontsize',17)

legend('X01 - X03 (AVG=0.655)','X03 - X04 (AVG=0.437)',...
    'X04 - X05 (AVG=-0.057)','X05 - X11 (AVG=1.284)','fontsize',10,...
    'location','north')
xlabel('Date','fontsize',15)
ylabel('Observed TED [W/m]','fontsize',15)
ylim([-3 11])
xlim([datetime(2022,6,15) datetime(2022,7,20)])

    % Asilomar - Sea
subplot(2,2,4)
plot(Tday,X.EWM_DD_OTEDSea{1},'bx','MarkerSize',10)
hold on; grid on;
plot(Tday,X.EWM_DD_OTEDSea{2},'mx','MarkerSize',10)
plot(TdayADCP,X.EWM_DD_OTEDSea{3},'gx','MarkerSize',10)
plot(TdayADCP,X.EWM_DD_OTEDSea{4},'rx','MarkerSize',10)
title('Asilomar: Sea Observed Energy Dissipation','fontsize',17)

legend('X01 - X03 (AVG=0.434)','X03 - X04 (AVG=0.425)',...
    'X04 - X05 (AVG=-0.545)','X05 - X11 (AVG=3.566)','fontsize',10,...
    'location','north')
xlabel('Date','fontsize',15)
ylabel('Observed TED [W/m]','fontsize',15)
ylim([-3 11])
xlim([datetime(2022,6,15) datetime(2022,7,20)])



%% Clean Up

clear i j N m0Sea m0Swell m1Sea m1Swell





















