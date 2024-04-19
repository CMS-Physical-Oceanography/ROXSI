%% plot_DailyNielsenTEDs.m
%
% Noah Clark        4/11/2024
%
%
% Purpose: To perform Nielsen TED analysis and determine the EWM per hour
%          in the sea and swell bands of these TED values. Then plot these
%          Nielsen TED daily averages. This is done for all spotter buoys
%          and ADCPs at China Rock and Asilomar.
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


%% Preliminaries

clc;clear;

% Load data for China Rock
load('WBvariables.mat','BSee','Btime','Bdepth')
% Load data for Asilomar
load('WBvariables.mat','XSee','Xtime','Xdepth')
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

a1 = 5.5;   a2 = -0.2;    a3 = -6.3;
kw = 0.5;

for i = 1:3
    BSee{i} = BSee{i}(2:end,:);
    XSee{i} = XSee{i}(2:end,:);
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



%% China Rock

        % DETERMINE ALL OF THE NIELSEN TEDS:
    % B01 to B03
% B01
[~,~,~,NTED1] = NC_NeilsonDiss(BSee{1},ff,Bdepth{1},kw,...
    a1,a2,a3,'r');
% B03
[~,~,~,NTED2] = NC_NeilsonDiss(BSee{2}(:,[1:397 505:end]),...
    ff,Bdepth{2}([1:397 505:end]),kw,a1,a2,a3,'r');
% Average
B.AllNielsen{1} = (NTED1 + NTED2)./2;
clear NTED1 NTED2

    % B03 to B05
% B03
[~,~,~,NTED1] = NC_NeilsonDiss(BSee{2},ff,Bdepth{2},kw,...
    a1,a2,a3,'r');
% B05
[~,~,~,NTED2] = NC_NeilsonDiss(BSee{3},ff,Bdepth{3},kw,...
    a1,a2,a3,'r');
% Average
B.AllNielsen{2} = (NTED1 + NTED2)./2;
clear NTED1 NTED2

    % B05 to B10
% B05
[~,~,~,NTED1] = NC_NeilsonDiss(BSee{3},ff,Bdepth{3},kw,...
    a1,a2,a3,'r');
% B10
[~,~,~,NTED2] = NC_NeilsonDiss(ADCP.B10.See,ff(1:35),ADCP.B10.depth,kw,...
    a1,a2,a3,'r');
% Average
B.AllNielsen{3} = (NTED1(1:35,217:817) + NTED2)./2;
clear NTED1 NTED2

    % B10 to B13
% B10
[~,~,~,NTED1] = NC_NeilsonDiss(ADCP.B10.See,ff(1:35),ADCP.B10.depth,kw,...
    a1,a2,a3,'r');
% B13
[~,~,~,NTED2] = NC_NeilsonDiss(ADCP.B13.See,ff(1:35),ADCP.B13.depth,kw,...
    a1,a2,a3,'r');
% Average
B.AllNielsen{4} = (NTED1 + NTED2)./2;
clear NTED1 NTED2



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
        m1Sea = trapz(ff(9:24),B.ASee{j}(9:24,i).*B.AllNielsen{j}(9:24,i));
        B.EWM_Nielsen_Sea{j}(i) = m1Sea./m0Sea;
            % Swell
        m0Swell = trapz(ff(4:9),B.ASee{j}(4:9,i),1);
        m1Swell = trapz(ff(4:9),B.ASee{j}(4:9,i).*B.AllNielsen{j}(4:9,i));
        B.EWM_Nielsen_Swell{j}(i) = m1Swell./m0Swell;
        
    end
end


        % AVERAGE PER DAY
% B01 to B03
% (**1st Section**)
for i = 1:15
    EWM_DD_Nielsen_Sea1(i) = mean(B.EWM_Nielsen_Sea{1}(i*24-23:i*24));
    EWM_DD_Nielsen_Swell1(i) = mean(B.EWM_Nielsen_Swell{1}(i*24-23:i*24));
end
% (**2nd Section**)
for i = 17.625:29.625
    EWM_DD_Nielsen_Sea2(i-16.625) = mean(B.EWM_Nielsen_Sea{1}(i*24-23:i*24));
    EWM_DD_Nielsen_Swell2(i-16.625) = mean(B.EWM_Nielsen_Swell{1}(i*24-23:i*24));
end

B.EWM_DD_NielsenSea{1} = cat(2,EWM_DD_Nielsen_Sea1,EWM_DD_Nielsen_Sea2);
B.EWM_DD_NielsenSwell{1} = cat(2,EWM_DD_Nielsen_Swell1,EWM_DD_Nielsen_Swell2);
clear EWM_DD_Nielsen_Sea1 EWM_DD_Nielsen_Sea2
clear EWM_DD_Nielsen_Swell1 EWM_DD_Nielsen_Swell2

% B03 to B05, B05 to B10, & B10 to B13
for j = 2:4
    if j == 2 % B03 to B05
        N = 34;
    else % B05 to B10 & B10 to B13
        N = 25;
    end
    
    for i = 1:N
        B.EWM_DD_NielsenSea{j}(i) = mean(B.EWM_Nielsen_Sea{j}(i*24-23:i*24));
        B.EWM_DD_NielsenSwell{j}(i) = mean(B.EWM_Nielsen_Swell{j}(i*24-23:i*24));
    end
end



%% Asilomar

        % DETERMINE ALL OF THE NIELSEN TEDS:
    % X01 to X03
% X01
[~,~,~,NTED1] = NC_NeilsonDiss(XSee{1},ff,Xdepth{1},kw,...
    a1,a2,a3,'r');
% X03
[~,~,~,NTED2] = NC_NeilsonDiss(XSee{2},ff,Xdepth{2},kw,...
    a1,a2,a3,'r');
% Average
X.AllNielsen{1} = (NTED1 + NTED2)./2;
clear NTED1 NTED2

    % X03 to X04
% X03
[~,~,~,NTED1] = NC_NeilsonDiss(XSee{2},ff,Xdepth{2},kw,...
    a1,a2,a3,'r');
% X04
[~,~,~,NTED2] = NC_NeilsonDiss(XSee{3},ff,Xdepth{3},kw,...
    a1,a2,a3,'r');
% Average
X.AllNielsen{2} = (NTED1 + NTED2)./2;
clear NTED1 NTED2

    % X04 to X05
% X04
[~,~,~,NTED1] = NC_NeilsonDiss(XSee{3},ff,Xdepth{3},kw,...
    a1,a2,a3,'r');
% X05
[~,~,~,NTED2] = NC_NeilsonDiss(ADCP.X05.See,ff(1:35),ADCP.X05.depth,kw,...
    a1,a2,a3,'r');
% Average
X.AllNielsen{3} = (NTED1(1:35,217:817) + NTED2)./2;
clear NTED1 NTED2

    % X05 to X11
% X05
[~,~,~,NTED1] = NC_NeilsonDiss(ADCP.X05.See,ff(1:35),ADCP.X05.depth,kw,...
    a1,a2,a3,'r');
% X11
[~,~,~,NTED2] = NC_NeilsonDiss(ADCP.X11.See,ff(1:35),ADCP.X11.depth,kw,...
    a1,a2,a3,'r');
% Average
X.AllNielsen{4} = (NTED1 + NTED2)./2;
clear NTED1 NTED2



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
        m1Sea = trapz(ff(9:24),X.ASee{j}(9:24,i).*X.AllNielsen{j}(9:24,i));
        X.EWM_Nielsen_Sea{j}(i) = m1Sea./m0Sea;
            % Swell
        m0Swell = trapz(ff(4:9),X.ASee{j}(4:9,i),1);
        m1Swell = trapz(ff(4:9),X.ASee{j}(4:9,i).*X.AllNielsen{j}(4:9,i));
        X.EWM_Nielsen_Swell{j}(i) = m1Swell./m0Swell;
        
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
        X.EWM_DD_NielsenSea{j}(i) = nanmean(X.EWM_Nielsen_Sea{j}(i*24-23:i*24));
        X.EWM_DD_NielsenSwell{j}(i) = nanmean(X.EWM_Nielsen_Swell{j}(i*24-23:i*24));
    end
end



%% Average Nielsen TED b/t Buoys

for i = 1:4
    % China Rock
    B.Avg_NTED_Sea{i} = nanmean(B.EWM_DD_NielsenSea{i});
    B.Avg_NTED_Swell{i} = nanmean(B.EWM_DD_NielsenSwell{i});
    % Asilomar
    X.Avg_NTED_Sea{i} = nanmean(X.EWM_DD_NielsenSea{i});
    X.Avg_NTED_Swell{i} = nanmean(X.EWM_DD_NielsenSwell{i});
end



%% Plotting

figure(1);clf;

    % China Rock - SWELL
subplot(2,2,1)
plot(TdayB01,B.EWM_DD_NielsenSwell{1},'bx','MarkerSize',10)
hold on; grid on;
plot(Tday,B.EWM_DD_NielsenSwell{2},'mx','MarkerSize',10)
plot(TdayADCP,B.EWM_DD_NielsenSwell{3},'gx','MarkerSize',10)
plot(TdayADCP,B.EWM_DD_NielsenSwell{4},'rx','MarkerSize',10)

title('China Rock: Swell Nielsen Energy Dissipation','fontsize',17)
legend('B01 - B03 (AVG=0.117)','B03 - B05 (AVG=0.192)',...
    'B05 - B10 (AVG=0.306)','B10 - B13 (AVG=0.727)','fontsize',10,...
    'location','north')
xlabel('Date','fontsize',15)
ylabel('Nielsen TED [W/m]','fontsize',15)
ylim([0 2.5])
xlim([datetime(2022,6,15) datetime(2022,7,20)])

    % China Rock - SEA
subplot(2,2,2)
plot(TdayB01,B.EWM_DD_NielsenSea{1},'bx','MarkerSize',10)
hold on; grid on;
plot(Tday,B.EWM_DD_NielsenSea{2},'mx','MarkerSize',10)
plot(TdayADCP,B.EWM_DD_NielsenSea{3},'gx','MarkerSize',10)
plot(TdayADCP,B.EWM_DD_NielsenSea{4},'rx','MarkerSize',10)

title('China Rock: Sea Nielsen Energy Dissipation','fontsize',17)
legend('B01 - B03 (AVG=0.111)','B03 - B05 (AVG=0.204)',...
    'B05 - B10 (AVG=0.451)','B10 - B13 (AVG=0.900)','fontsize',10,...
    'location','north')
xlabel('Date','fontsize',15)
ylabel('Nielsen TED [W/m]','fontsize',15)
ylim([0 2.5])
xlim([datetime(2022,6,15) datetime(2022,7,20)])

    % Asilomar - SWELL
subplot(2,2,3)
plot(Tday,X.EWM_DD_NielsenSwell{1},'bx','MarkerSize',10)
hold on; grid on;
plot(Tday,X.EWM_DD_NielsenSwell{2},'mx','MarkerSize',10)
plot(TdayADCP,X.EWM_DD_NielsenSwell{3},'gx','MarkerSize',10)
plot(TdayADCP,X.EWM_DD_NielsenSwell{4},'rx','MarkerSize',10)

title('Asilomar: Swell Nielsen Energy Dissipation','fontsize',17)
legend('X01 - X03 (AVG=0.177)','X03 - X04 (AVG=0.262)',...
    'X04 - X05 (AVG=0.206)','X05 - X11 (AVG=0.279)','fontsize',10,...
    'location','north')
xlabel('Date','fontsize',15)
ylabel('Nielsen TED [W/m]','fontsize',15)
ylim([0 1.65])
xlim([datetime(2022,6,15) datetime(2022,7,20)])

    % Asilomar - SEA
subplot(2,2,4)
plot(Tday,X.EWM_DD_NielsenSea{1},'bx','MarkerSize',10)
hold on; grid on;
plot(Tday,X.EWM_DD_NielsenSea{2},'mx','MarkerSize',10)
plot(TdayADCP,X.EWM_DD_NielsenSea{3},'gx','MarkerSize',10)
plot(TdayADCP,X.EWM_DD_NielsenSea{4},'rx','MarkerSize',10)

title('Asilomar: Sea Nielsen Energy Dissipation','fontsize',17)
legend('X01 - X03 (AVG=0.241)','X03 - X04 (AVG=0.431)',...
    'X04 - X05 (AVG=0.601)','X05 - X11 (AVG=0.513)','fontsize',10,...
    'location','north')
xlabel('Date','fontsize',15)
ylabel('Nielsen TED [W/m]','fontsize',15)
ylim([0 1.65])
xlim([datetime(2022,6,15) datetime(2022,7,20)])



%% Clean Up

clear i j N m0Sea m0Swell m1Sea m1Swell





