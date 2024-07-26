%% plot_ErrBarAvgSpecs.m
%
% Noah Clark        1/11/24
%
%
% Purpose: Determine 3 time periods of 1 day and 15 hours where there is
%          only a peak in energy in the sea band, swell band, and where 
%          there is a peak at both. Plot the energy, wave direction, and
%          directional spread during these times with a std enelope over
%          each. 
%          Then compare the obs TED and Nielsen TED for these 3 time
%          periods. 
%
% ------------------------------------------------------------------------
%
% Main Variables:
%                - SeaP: structure containing all of the data used to make
%                        the plots for the selected sea peak time period
%                - SwellP: structure containing all of the data used to
%                          make the plots for the selected swell peak time
%                          period 
%                - BothP : structure containing all of the data used to
%                          make the plots for the selected time period
%                          where there is a peak at both the sea and swell
%                          frequency ranges
%
%
% Notes:
%           - Swell Range: 0.049 - 0.088 Hz     (ff(5)-ff(9))
%           - Sea Range: 0.088 - 0.137 Hz       (ff(9)-ff(14))
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


%% General Preliminaries
% The general preliminaries for the original plots of the 3 selected time
% periods at each of the 3 spotter buoys (9 plots)

clc;clear;

addpath('../Summer (0516-0728)/Start Data')

% Load Normal Wave Direction: 
load('WBvariables.mat','BNormWaveDir','Bmeanspec','BEMEMmean','Btime','BSee')

    % Define the difference between Sea and Swell:
Def_Freq = 0.088; % >: Sea       % <: Swell

    % Frequency Vector: 
load('roxsi_spotter_L2_B05_1153_reduced.mat') %just load any to extract freq
ff = spotterL2.frequency;


    % Define the three 24 hour time periods (Sea, Peak, & Both):
% Sea Period: [July 6 @ 12:00 - July 7 @ 12:00]   (1 day)
SeaP.ind1 = 505;
SeaP.ind2 = 529;
SeaP.N = SeaP.ind1 - SeaP.ind2 + 1;
SeaP.TimeStart = spotterL2.dtime(SeaP.ind1); % July 6 @ 12:00
SeaP.TimeEnd = spotterL2.dtime(SeaP.ind2); % July 7 @ 12:00
% Swell Period: [June 26 @ 04:00 - June 27 @ 04:00]      (1 day)
SwellP.ind1 = 257;
SwellP.ind2 = 281;
SwellP.N = SwellP.ind1 - SwellP.ind2 + 1;
SwellP.TimeStart = spotterL2.dtime(SwellP.ind1); % June 26 @ 04:00
SwellP.TimeEnd = spotterL2.dtime(SwellP.ind2); % June 27 @ 04:00
% Both Period: [Jun 30 @ 20:00 - July 01 @ 20:00]      (1 day)
BothP.ind1 = 365;
BothP.ind2 = 389;
BothP.N = BothP.ind1 - BothP.ind2 + 1;
BothP.TimeStart = spotterL2.dtime(BothP.ind1); % June 30 @ 20:00
BothP.TimeEnd = spotterL2.dtime(BothP.ind2); % July 01 @ 20:00


%% Spectrogram of B05 with Marked Periods

% The y locations
yL = [0 0.35];

% Sea Period x locations
Seax1 = [datetime(string(SeaP.TimeStart),'TimeZone','America/Los_Angeles')...
        , datetime(string(SeaP.TimeStart),'TimeZone','America/Los_Angeles')];
Seax2 = [datetime(string(SeaP.TimeEnd),'TimeZone','America/Los_Angeles')...
        , datetime(string(SeaP.TimeEnd),'TimeZone','America/Los_Angeles')];
% Swell Period x locations
Swellx1 = [datetime(string(SwellP.TimeStart),'TimeZone','America/Los_Angeles')...
        , datetime(string(SwellP.TimeStart),'TimeZone','America/Los_Angeles')];
Swellx2 = [datetime(string(SwellP.TimeEnd),'TimeZone','America/Los_Angeles')...
        , datetime(string(SwellP.TimeEnd),'TimeZone','America/Los_Angeles')];
% Both Period x locations
Bothx1 = [datetime(string(BothP.TimeStart),'TimeZone','America/Los_Angeles')...
        , datetime(string(BothP.TimeStart),'TimeZone','America/Los_Angeles')];
Bothx2 = [datetime(string(BothP.TimeEnd),'TimeZone','America/Los_Angeles')...
        , datetime(string(BothP.TimeEnd),'TimeZone','America/Los_Angeles')];

figure(20);clf;
pcolor(Btime{3},ff,BSee{3})
hold on
plot(Swellx1,yL,'r','LineWidth',1.5)
plot(Swellx2,yL,'r','LineWidth',1.5)
plot(Bothx1,yL,'g','LineWidth',1.5)
plot(Bothx2,yL,'g','LineWidth',1.5)
plot(Seax1,yL,'m','LineWidth',1.5)
plot(Seax2,yL,'m','LineWidth',1.5)

ylim([0 0.35])
title('B05 Spectrogram','fontsize',18)
shading flat
cb = colorbar;
caxis([0 3])
xlabel('Date','fontsize',16)
ylabel(cb,'Wave Energy (m^2/Hz)','fontsize',16)
ylabel('Frequency (Hz)','fontsize',16)
legend('','Swell Period','','Both Period','','Sea Period','fontsize',14)


%% B01 Preliminaries

    % Fully Time Averaged Variables (For B01):
load('roxsi_spotter_L2_B01_1158_reduced.mat') % first part for B01
DS1 = spotterL2.EMEM.meandirspread_f;
load('roxsi_spotter_L2_B01_1150_reduced.mat') % second part for B01
DS2 = spotterL2.EMEM.meandirspread_f;
cDS = cat(2,DS1,DS2);

B01.FTA_See = Bmeanspec{1};
B01.FTA_Dir = BEMEMmean{1};
B01.FTA_DirSpread = meanangle(cDS,2);

clear DS1 DS2 cDS


%% B01: Sea Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy only peaks in the sea band

% Load in 2nd part of B01:
load('roxsi_spotter_L2_B01_1150_reduced.mat')

    % Assign Variables (for the day and a half):
% Energy
SeaP.B01.mSee = mean(spotterL2.See(:,SeaP.ind1-504:SeaP.ind2-504),2);
SeaP.B01.std_mSee = std(spotterL2.See(:,SeaP.ind1-504:SeaP.ind2-504),0,2);
SeaP.B01.MMSee = mean(SeaP.B01.mSee(9:14));
% Direction
SeaP.B01.mDir = meanangle(spotterL2.EMEM.meandir_f(:,SeaP.ind1-504:SeaP.ind2-504),2) + 360;
SeaP.B01.std_mDir = std(spotterL2.EMEM.meandir_f(:,SeaP.ind1-504:SeaP.ind2-504),0,2);
SeaP.B01.MMDir = mean(SeaP.B01.mDir(9:14));
% Directional Spread
SeaP.B01.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,SeaP.ind1-504:SeaP.ind2-504),2);
SeaP.B01.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,SeaP.ind1-504:SeaP.ind2-504),0,2);
SeaP.B01.MMDirSpread = mean(SeaP.B01.mDirSpread(9:14));

    % Plot:
figure(1);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,SeaP.B01.mSee,SeaP.B01.std_mSee,'k2')
hold on
plot(ff,SeaP.B01.mSee,'r','LineWidth',1.5)
plot(ff,B01.FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.83,0.1,0.05],'string',sprintf('%0.3f',SeaP.B01.MMSee));
grid on
title('B01 Sea Peak Period: Mean Energy (Jul 6 @ 12:00 - Jul 7 @ 12:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.5 1.3])
xlim([0 ff(37)])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,SeaP.B01.mDir,SeaP.B01.std_mDir,'k2')
hold on
plot(ff,SeaP.B01.mDir,'b','LineWidth',1.5)
plot(ff,B01.FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
an = annotation('textbox',[0.26,0.53,0.1,0.05],'string',sprintf('%0.1f',SeaP.B01.MMDir));
grid on
title('B01 Sea Peak Period: Mean Direction (Jul 6 @ 12:00 - Jul 7 @ 12:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([225 325])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','east')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,SeaP.B01.mDirSpread,SeaP.B01.std_mDirSpread,'k2')
hold on
plot(ff,SeaP.B01.mDirSpread,'g','LineWidth',1.5)
plot(ff,B01.FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.2,0.1,0.05],'string',sprintf('%0.1f',SeaP.B01.MMDirSpread));
grid on
title('B01 Sea Peak Period: Mean Directional Spread (Jul 6 @ 12:00 - Jul 7 @ 12:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% B01: Swell Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy only peaks in the swell band

% Load in 1st part of B01:
load('roxsi_spotter_L2_B01_1158_reduced.mat')

    % Assign Variables (for the day and a half):
% Energy
SwellP.B01.mSee = mean(spotterL2.See(:,SwellP.ind1:SwellP.ind2),2);
SwellP.B01.std_mSee = std(spotterL2.See(:,SwellP.ind1:SwellP.ind2),0,2);
SwellP.B01.MMSee = mean(SwellP.B01.mSee(5:9));
% Direction
SwellP.B01.mDir = meanangle(spotterL2.EMEM.meandir_f(:,SwellP.ind1:SwellP.ind2),2) + 360;
SwellP.B01.std_mDir = std(spotterL2.EMEM.meandir_f(:,SwellP.ind1:SwellP.ind2),0,2);
SwellP.B01.MMDir = mean(SwellP.B01.mDir(5:9));
% Directional Spread
SwellP.B01.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,SwellP.ind1:SwellP.ind2),2);
SwellP.B01.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,SwellP.ind1:SwellP.ind2),0,2);
SwellP.B01.MMDirSpread = mean(SwellP.B01.mDirSpread(5:9));

    % Plot:
figure(2);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,SwellP.B01.mSee,SwellP.B01.std_mSee,'k2')
hold on
plot(ff,SwellP.B01.mSee,'r','LineWidth',1.5)
plot(ff,B01.FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.83,0.1,0.05],'string',sprintf('%0.3f',SwellP.B01.MMSee));
grid on
title('B01 Swell Peak Period: Mean Energy (Jun 26 @ 04:00 - Jun 27 @ 04:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.25 1.2])
xlim([0 ff(37)])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,SwellP.B01.mDir,SwellP.B01.std_mDir,'k2')
hold on
plot(ff,SwellP.B01.mDir,'b','LineWidth',1.5)
plot(ff,B01.FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
an = annotation('textbox',[0.26,0.53,0.1,0.05],'string',sprintf('%0.1f',SwellP.B01.MMDir));
grid on
title('B01 Swell Peak Period: Mean Direction (Jun 26 @ 04:00 - Jun 27 @ 04:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([175 325])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','east')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,SwellP.B01.mDirSpread,SwellP.B01.std_mDirSpread,'k2')
hold on
plot(ff,SwellP.B01.mDirSpread,'g','LineWidth',1.5)
plot(ff,B01.FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.2,0.1,0.05],'string',sprintf('%0.1f',SwellP.B01.MMDirSpread));
grid on
title('B01 Swell Peak Period: Mean Directional Spread (Jun 26 @ 04:00 - Jun 27 @ 04:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% B01: Both Swell and Sea Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy peaks in both the swell and the sea bands

    % Assign Variables (for the day and a half):
% Energy
BothP.B01.mSee = mean(spotterL2.See(:,BothP.ind1:BothP.ind2),2);
BothP.B01.std_mSee = std(spotterL2.See(:,BothP.ind1:BothP.ind2),0,2);
BothP.B01.MMSee1 = mean(BothP.B01.mSee(5:9));
BothP.B01.MMSee2 = mean(BothP.B01.mSee(9:14));
% Direction
BothP.B01.mDir = meanangle(spotterL2.EMEM.meandir_f(:,BothP.ind1:BothP.ind2),2) + 360;
BothP.B01.std_mDir = std(spotterL2.EMEM.meandir_f(:,BothP.ind1:BothP.ind2),0,2);
BothP.B01.MMDir1 = mean(BothP.B01.mDir(5:9));
BothP.B01.MMDir2 = mean(BothP.B01.mDir(9:14));
% Directional Spread
BothP.B01.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,BothP.ind1:BothP.ind2),2);
BothP.B01.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,BothP.ind1:BothP.ind2),0,2);
BothP.B01.MMDirSpread1 = mean(BothP.B01.mDirSpread(5:9));
BothP.B01.MMDirSpread2 = mean(BothP.B01.mDirSpread(9:14));

    % Plot:
figure(3);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,BothP.B01.mSee,BothP.B01.std_mSee,'k2')
hold on
plot(ff,BothP.B01.mSee,'r','LineWidth',1.5)
plot(ff,B01.FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an1 = annotation('textbox',[0.26,0.83,0.1,0.05],'string',sprintf('%0.3f',BothP.B01.MMSee1));
an2 = annotation('textbox',[0.44,0.83,0.1,0.05],'string',sprintf('%0.3f',BothP.B01.MMSee2));
grid on
title('B01 Both Peaks Period: Mean Energy (Jun 30 @ 20:00 - Jul 01 @ 20:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.5 2])
xlim([0 ff(37)])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,BothP.B01.mDir,BothP.B01.std_mDir,'k2')
hold on
plot(ff,BothP.B01.mDir,'b','LineWidth',1.5)
plot(ff,B01.FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
an1 = annotation('textbox',[0.26,0.53,0.1,0.05],'string',sprintf('%0.1f',BothP.B01.MMDir1));
an2 = annotation('textbox',[0.44,0.47,0.1,0.05],'string',sprintf('%0.1f',BothP.B01.MMDir2));
grid on
title('B01 Both Peaks Period: Mean Direction (Jun 30 @ 20:00 - Jul 01 @ 20:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([225 325])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','east')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,BothP.B01.mDirSpread,BothP.B01.std_mDirSpread,'k2')
hold on
plot(ff,BothP.B01.mDirSpread,'g','LineWidth',1.5)
plot(ff,B01.FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an1 = annotation('textbox',[0.26,0.22,0.1,0.05],'string',sprintf('%0.1f',BothP.B01.MMDirSpread1));
an2 = annotation('textbox',[0.44,0.22,0.1,0.05],'string',sprintf('%0.1f',BothP.B01.MMDirSpread2));
grid on
title('B01 Both Peaks Period: Mean Directional Spread (Jun 30 @ 20:00 - Jul 01 @ 20:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell','location','east')


%% B03 Preliminaries

    % Fully Time Averaged Variables (For B03):
load('roxsi_spotter_L2_B03_1152_reduced.mat')
B03.FTA_See = mean(spotterL2.See,2);
B03.FTA_Dir = meanangle(spotterL2.EMEM.meandir_f,2) + 360;
B03.FTA_DirSpread = meanangle(spotterL2.EMEM.meandirspread_f,2);


%% B03: Sea Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy only peaks in the sea band

    % Assign Variables (for the day and a half):
% Energy
SeaP.B03.mSee = mean(spotterL2.See(:,SeaP.ind1:SeaP.ind2),2);
SeaP.B03.std_mSee = std(spotterL2.See(:,SeaP.ind1:SeaP.ind2),0,2);
SeaP.B03.MMSee = mean(SeaP.B03.mSee(9:14));
% Direction
SeaP.B03.mDir = meanangle(spotterL2.EMEM.meandir_f(:,SeaP.ind1:SeaP.ind2),2) + 360;
SeaP.B03.std_mDir = std(spotterL2.EMEM.meandir_f(:,SeaP.ind1:SeaP.ind2),0,2);
SeaP.B03.MMDir = mean(SeaP.B03.mDir(9:14));
% Directional Spread
SeaP.B03.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,SeaP.ind1:SeaP.ind2),2);
SeaP.B03.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,SeaP.ind1:SeaP.ind2),0,2);
SeaP.B03.MMDirSpread = mean(SeaP.B03.mDirSpread(9:14));

    % Plot:
figure(4);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,SeaP.B03.mSee,SeaP.B03.std_mSee,'k2')
hold on
plot(ff,SeaP.B03.mSee,'r','LineWidth',1.5)
plot(ff,B03.FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.83,0.1,0.05],'string',sprintf('%0.3f',SeaP.B03.MMSee));
grid on
title('B03 Sea Peak Period: Mean Energy (Jul 6 @ 12:00 - Jul 7 @ 12:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.5 1.3])
xlim([0 ff(37)])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,SeaP.B03.mDir,SeaP.B03.std_mDir,'k2')
hold on
plot(ff,SeaP.B03.mDir,'b','LineWidth',1.5)
plot(ff,B03.FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
an = annotation('textbox',[0.26,0.53,0.1,0.05],'string',sprintf('%0.1f',SeaP.B03.MMDir));
grid on
title('B03 Sea Peak Period: Mean Direction (Jul 6 @ 12:00 - Jul 7 @ 12:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([225 325])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','east')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,SeaP.B03.mDirSpread,SeaP.B03.std_mDirSpread,'k2')
hold on
plot(ff,SeaP.B03.mDirSpread,'g','LineWidth',1.5)
plot(ff,B03.FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.2,0.1,0.05],'string',sprintf('%0.1f',SeaP.B03.MMDirSpread));
grid on
title('B03 Sea Peak Period: Mean Directional Spread (Jul 6 @ 12:00 - Jul 7 @ 12:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% B03: Swell Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy only peaks in the swell band

    % Assign Variables (for the day and a half):
% Energy
SwellP.B03.mSee = mean(spotterL2.See(:,SwellP.ind1:SwellP.ind2),2);
SwellP.B03.std_mSee = std(spotterL2.See(:,SwellP.ind1:SwellP.ind2),0,2);
SwellP.B03.MMSee = mean(SwellP.B03.mSee(5:9));
% Direction
SwellP.B03.mDir = meanangle(spotterL2.EMEM.meandir_f(:,SwellP.ind1:SwellP.ind2),2) + 360;
SwellP.B03.std_mDir = std(spotterL2.EMEM.meandir_f(:,SwellP.ind1:SwellP.ind2),0,2);
SwellP.B03.MMDir = mean(SwellP.B03.mDir(5:9));
% Directional Spread
SwellP.B03.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,SwellP.ind1:SwellP.ind2),2);
SwellP.B03.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,SwellP.ind1:SwellP.ind2),0,2);
SwellP.B03.MMDirSpread = mean(SwellP.B03.mDirSpread(5:9));

    % Plot:
figure(5);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,SwellP.B03.mSee,SwellP.B03.std_mSee,'k2')
hold on
plot(ff,SwellP.B03.mSee,'r','LineWidth',1.5)
plot(ff,B03.FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.83,0.1,0.05],'string',sprintf('%0.3f',SwellP.B03.MMSee));
grid on
title('B03 Swell Peak Period: Mean Energy (Jun 26 @ 04:00 - Jun 27 @ 04:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.25 1.2])
xlim([0 ff(37)])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,SwellP.B03.mDir,SwellP.B03.std_mDir,'k2')
hold on
plot(ff,SwellP.B03.mDir,'b','LineWidth',1.5)
plot(ff,B03.FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
an = annotation('textbox',[0.26,0.53,0.1,0.05],'string',sprintf('%0.1f',SwellP.B03.MMDir));
grid on
title('B03 Swell Peak Period: Mean Direction (Jun 26 @ 04:00 - Jun 27 @ 04:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([175 325])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','east')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,SwellP.B03.mDirSpread,SwellP.B03.std_mDirSpread,'k2')
hold on
plot(ff,SwellP.B03.mDirSpread,'g','LineWidth',1.5)
plot(ff,B03.FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.2,0.1,0.05],'string',sprintf('%0.1f',SwellP.B03.MMDirSpread));
grid on
title('B03 Swell Peak Period: Mean Directional Spread (Jun 26 @ 04:00 - Jun 27 @ 04:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% B03: Both Swell and Sea Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy peaks in both the swell and the sea bands

    % Assign Variables (for the day and a half):
% Energy
BothP.B03.mSee = mean(spotterL2.See(:,BothP.ind1:BothP.ind2),2);
BothP.B03.std_mSee = std(spotterL2.See(:,BothP.ind1:BothP.ind2),0,2);
BothP.B03.MMSee1 = mean(BothP.B03.mSee(5:9));
BothP.B03.MMSee2 = mean(BothP.B03.mSee(9:14));
% Direction
BothP.B03.mDir = meanangle(spotterL2.EMEM.meandir_f(:,BothP.ind1:BothP.ind2),2) + 360;
BothP.B03.std_mDir = std(spotterL2.EMEM.meandir_f(:,BothP.ind1:BothP.ind2),0,2);
BothP.B03.MMDir1 = mean(BothP.B03.mDir(5:9));
BothP.B03.MMDir2 = mean(BothP.B03.mDir(9:14));
% Directional Spread
BothP.B03.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,BothP.ind1:BothP.ind2),2);
BothP.B03.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,BothP.ind1:BothP.ind2),0,2);
BothP.B03.MMDirSpread1 = mean(BothP.B03.mDirSpread(5:9));
BothP.B03.MMDirSpread2 = mean(BothP.B03.mDirSpread(9:14));

    % Plot:
figure(6);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,BothP.B03.mSee,BothP.B03.std_mSee,'k2')
hold on
plot(ff,BothP.B03.mSee,'r','LineWidth',1.5)
plot(ff,B03.FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an1 = annotation('textbox',[0.26,0.83,0.1,0.05],'string',sprintf('%0.3f',BothP.B03.MMSee1));
an2 = annotation('textbox',[0.44,0.83,0.1,0.05],'string',sprintf('%0.3f',BothP.B03.MMSee2));
grid on
title('B03 Both Peaks Period: Mean Energy (Jun 30 @ 20:00 - Jul 01 @ 20:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.5 2])
xlim([0 ff(37)])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,BothP.B03.mDir,BothP.B03.std_mDir,'k2')
hold on
plot(ff,BothP.B03.mDir,'b','LineWidth',1.5)
plot(ff,B03.FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
an1 = annotation('textbox',[0.26,0.53,0.1,0.05],'string',sprintf('%0.1f',BothP.B03.MMDir1));
an2 = annotation('textbox',[0.44,0.47,0.1,0.05],'string',sprintf('%0.1f',BothP.B03.MMDir2));
grid on
title('B03 Both Peaks Period: Mean Direction (Jun 30 @ 20:00 - Jul 01 @ 20:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([225 325])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','east')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,BothP.B03.mDirSpread,BothP.B03.std_mDirSpread,'k2')
hold on
plot(ff,BothP.B03.mDirSpread,'g','LineWidth',1.5)
plot(ff,B03.FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an1 = annotation('textbox',[0.26,0.22,0.1,0.05],'string',sprintf('%0.1f',BothP.B03.MMDirSpread1));
an2 = annotation('textbox',[0.44,0.22,0.1,0.05],'string',sprintf('%0.1f',BothP.B03.MMDirSpread2));
grid on
title('B03 Both Peaks Period: Mean Directional Spread (Jun 30 @ 20:00 - Jul 01 @ 20:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% B05 Preliminaries

    % Fully Time Averaged Variables (For B05):
load('roxsi_spotter_L2_B05_1153_reduced.mat')
B05.FTA_See = mean(spotterL2.See,2);
B05.FTA_Dir = meanangle(spotterL2.EMEM.meandir_f,2) + 360;
B05.FTA_DirSpread = meanangle(spotterL2.EMEM.meandirspread_f,2);


%% B05: Sea Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy only peaks in the sea band

    % Assign Variables (for the day and a half):
% Energy
SeaP.B05.mSee = mean(spotterL2.See(:,SeaP.ind1:SeaP.ind2),2);
SeaP.B05.std_mSee = std(spotterL2.See(:,SeaP.ind1:SeaP.ind2),0,2);
SeaP.B05.MMSee = mean(SeaP.B05.mSee(9:14));
% Direction
SeaP.B05.mDir = meanangle(spotterL2.EMEM.meandir_f(:,SeaP.ind1:SeaP.ind2),2) + 360;
SeaP.B05.std_mDir = std(spotterL2.EMEM.meandir_f(:,SeaP.ind1:SeaP.ind2),0,2);
SeaP.B05.MMDir = mean(SeaP.B05.mDir(9:14));
% Directional Spread
SeaP.B05.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,SeaP.ind1:SeaP.ind2),2);
SeaP.B05.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,SeaP.ind1:SeaP.ind2),0,2);
SeaP.B05.MMDirSpread = mean(SeaP.B05.mDirSpread(9:14));

    % Plot:
figure(7);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,SeaP.B05.mSee,SeaP.B05.std_mSee,'k2')
hold on
plot(ff,SeaP.B05.mSee,'r','LineWidth',1.5)
plot(ff,B05.FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.83,0.1,0.05],'string',sprintf('%0.3f',SeaP.B05.MMSee));
grid on
title('B05 Sea Peak Period: Mean Energy (Jul 6 @ 12:00 - Jul 7 @ 12:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.5 1.3])
xlim([0 ff(37)])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,SeaP.B05.mDir,SeaP.B05.std_mDir,'k2')
hold on
plot(ff,SeaP.B05.mDir,'b','LineWidth',1.5)
plot(ff,B05.FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
an = annotation('textbox',[0.26,0.53,0.1,0.05],'string',sprintf('%0.1f',SeaP.B05.MMDir));
grid on
title('B05 Sea Peak Period: Mean Direction (Jul 6 @ 12:00 - Jul 7 @ 12:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([225 325])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','east')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,SeaP.B05.mDirSpread,SeaP.B05.std_mDirSpread,'k2')
hold on
plot(ff,SeaP.B05.mDirSpread,'g','LineWidth',1.5)
plot(ff,B05.FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.2,0.1,0.05],'string',sprintf('%0.1f',SeaP.B05.MMDirSpread));
grid on
title('B05 Sea Peak Period: Mean Directional Spread (Jul 6 @ 12:00 - Jul 7 @ 12:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% B05: Swell Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy only peaks in the swell band

    % Assign Variables (for the day and a half):
% Energy
SwellP.B05.mSee = mean(spotterL2.See(:,SwellP.ind1:SwellP.ind2),2);
SwellP.B05.std_mSee = std(spotterL2.See(:,SwellP.ind1:SwellP.ind2),0,2);
SwellP.B05.MMSee = mean(SwellP.B05.mSee(5:9));
% Direction
SwellP.B05.mDir = meanangle(spotterL2.EMEM.meandir_f(:,SwellP.ind1:SwellP.ind2),2) + 360;
SwellP.B05.std_mDir = std(spotterL2.EMEM.meandir_f(:,SwellP.ind1:SwellP.ind2),0,2);
SwellP.B05.MMDir = mean(SwellP.B05.mDir(5:9));
% Directional Spread
SwellP.B05.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,SwellP.ind1:SwellP.ind2),2);
SwellP.B05.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,SwellP.ind1:SwellP.ind2),0,2);
SwellP.B05.MMDirSpread = mean(SwellP.B05.mDirSpread(5:9));

    % Plot:
figure(8);clf; 
%Energy
subplot(3,1,1)
bgpatch2(ff,SwellP.B05.mSee,SwellP.B05.std_mSee,'k2')
hold on
plot(ff,SwellP.B05.mSee,'r','LineWidth',1.5)
plot(ff,B05.FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.83,0.1,0.05],'string',sprintf('%0.3f',SwellP.B05.MMSee));
grid on
title('B05 Swell Peak Period: Mean Energy (Jun 26 @ 04:00 - Jun 27 @ 04:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.25 1.2])
xlim([0 ff(37)])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,SwellP.B05.mDir,SwellP.B05.std_mDir,'k2')
hold on
plot(ff,SwellP.B05.mDir,'b','LineWidth',1.5)
plot(ff,B05.FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
an = annotation('textbox',[0.26,0.53,0.1,0.05],'string',sprintf('%0.1f',SwellP.B05.MMDir));
grid on
title('B05 Swell Peak Period: Mean Direction (Jun 26 @ 04:00 - Jun 27 @ 04:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([175 325])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','east')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,SwellP.B05.mDirSpread,SwellP.B05.std_mDirSpread,'k2')
hold on
plot(ff,SwellP.B05.mDirSpread,'g','LineWidth',1.5)
plot(ff,B05.FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an = annotation('textbox',[0.44,0.2,0.1,0.05],'string',sprintf('%0.1f',SwellP.B05.MMDirSpread));
grid on
title('B05 Swell Peak Period: Mean Directional Spread (Jun 26 @ 04:00 - Jun 27 @ 04:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% B05: Both Swell and Sea Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy peaks in both the swell and the sea bands

    % Assign Variables (for the day and a half):
% Energy
BothP.B05.mSee = mean(spotterL2.See(:,BothP.ind1:BothP.ind2),2);
BothP.B05.std_mSee = std(spotterL2.See(:,BothP.ind1:BothP.ind2),0,2);
BothP.B05.MMSee1 = mean(BothP.B05.mSee(5:9));
BothP.B05.MMSee2 = mean(BothP.B05.mSee(9:14));
% Direction
BothP.B05.mDir = meanangle(spotterL2.EMEM.meandir_f(:,BothP.ind1:BothP.ind2),2) + 360;
BothP.B05.std_mDir = std(spotterL2.EMEM.meandir_f(:,BothP.ind1:BothP.ind2),0,2);
BothP.B05.MMDir1 = mean(BothP.B05.mDir(5:9));
BothP.B05.MMDir2 = mean(BothP.B05.mDir(9:14));
% Directional Spread
BothP.B05.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,BothP.ind1:BothP.ind2),2);
BothP.B05.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,BothP.ind1:BothP.ind2),0,2);
BothP.B05.MMDirSpread1 = mean(BothP.B05.mDirSpread(5:9));
BothP.B05.MMDirSpread2 = mean(BothP.B05.mDirSpread(9:14));

    % Plot:
figure(9);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,BothP.B05.mSee,BothP.B05.std_mSee,'k2')
hold on
plot(ff,BothP.B05.mSee,'r','LineWidth',1.5)
plot(ff,B05.FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an1 = annotation('textbox',[0.26,0.83,0.1,0.05],'string',sprintf('%0.3f',BothP.B05.MMSee1));
an2 = annotation('textbox',[0.44,0.83,0.1,0.05],'string',sprintf('%0.3f',BothP.B05.MMSee2));
grid on
title('B05 Both Peaks Period: Mean Energy (Jun 30 @ 20:00 - Jul 01 @ 20:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.5 2])
xlim([0 ff(37)])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,BothP.B05.mDir,BothP.B05.std_mDir,'k2')
hold on
plot(ff,BothP.B05.mDir,'b','LineWidth',1.5)
plot(ff,B05.FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
an1 = annotation('textbox',[0.26,0.53,0.1,0.05],'string',sprintf('%0.1f',BothP.B05.MMDir1));
an2 = annotation('textbox',[0.44,0.47,0.1,0.05],'string',sprintf('%0.1f',BothP.B05.MMDir2));
grid on
title('B05 Both Peaks Period: Mean Direction (Jun 30 @ 20:00 - Jul 01 @ 20:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([225 325])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','east')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,BothP.B05.mDirSpread,BothP.B05.std_mDirSpread,'k2')
hold on
plot(ff,BothP.B05.mDirSpread,'g','LineWidth',1.5)
plot(ff,B05.FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
an1 = annotation('textbox',[0.26,0.22,0.1,0.05],'string',sprintf('%0.1f',BothP.B05.MMDirSpread1));
an2 = annotation('textbox',[0.44,0.22,0.1,0.05],'string',sprintf('%0.1f',BothP.B05.MMDirSpread2));
grid on
title('B05 Both Peaks Period: Mean Directional Spread (Jun 30 @ 20:00 - Jul 01 @ 20:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries (for TED and Nielsen TED)

% Have to change ff from starting at 0 to starting at 0.0098 or else the
%  Nielsen function will return all NaNs:
ff = 0.0098:0.0098:(0.0098*129);
ffa = 0.0028:0.0028:(0.0028*126); %freq vector for original ADCP data


% Add the path for the functions
addpath('../Summer (0516-0728)/Main Functions')
addpath('../Fall (0825)')
addpath('../Fall (0825)/ADCP Data')

% Load data (1st set) for B01
load('roxsi_spotter_L2_B01_1158_reduced.mat')
See_B01_1st = spotterL2.See;
Direc_B01_1st = spotterL2.EMEM.meandir_f;
% Load data (2nd set) for B01 (NOT USED)
load('roxsi_spotter_L2_B01_1150_reduced.mat')
See_B01_2nd = spotterL2.See;
Direc_B01_2nd = spotterL2.EMEM.meandir_f;
% Load data for B03
load('roxsi_spotter_L2_B03_1152_reduced.mat')
See_B03 = spotterL2.See;
Direc_B03 = spotterL2.EMEM.meandir_f;
% Load data for B05
load('roxsi_spotter_L2_B05_1153_reduced.mat')
See_B05 = spotterL2.See;
Direc_B05 = spotterL2.EMEM.meandir_f;
% Load data for B10
load('roxsi_signature_L2_B10_See.mat')
See_B10 = A.See;
See_B10 = interp1(ffa,See_B10,ff);
See_B10 = See_B10(1:35,:);
Direc_B10 = Direc_B05(1:35,:);
Depth_B10 = A.depth;
% Load data for B13
load('roxsi_signature_L2_B13_See.mat')
See_B13 = A.See;
See_B13 = interp1(ffa,See_B13,ff);
See_B13 = See_B13(1:35,:);
Direc_B13 = Direc_B05(1:35,:);
Depth_B13 = A.depth;


% Load data needed for functions
load('WBvariables.mat','Butm','Bdepth')
load('SM&ADCP_All.mat')

% More needed for functions
kw = 0.5;   a1 = 5.5;   a2 = -0.2;    a3 = -6.3;


%% TEDs for Sea Peak Time Period

i1 = SeaP.ind1; 
i2 = SeaP.ind2;
for i = i1:i2
    % b/t B01 and B03:
    SeaP.ObsTED{1}(:,i+1-i1) = NC_ObsDiss(See_B01_2nd(:,i-504),See_B03(:,i),...
        Direc_B01_2nd(:,i-504),Direc_B03(:,i),ff,Butm{1},Butm{2},Bdepth{1}(i),...
        Bdepth{2}(i));
    [~,~,~,SeaP.NTED1{1}(:,i+1-i1)] = NC_NeilsonDiss(See_B01_2nd(:,i-504),ff,...
        Bdepth{1}(i),kw,a1,a2,a3,'r');
    [~,~,~,SeaP.NTED2{1}(:,i+1-i1)] = NC_NeilsonDiss(See_B03(:,i),ff,...
        Bdepth{2}(i),kw,a1,a2,a3,'r');
    
    % b/t B03 and B05:
    SeaP.ObsTED{2}(:,i+1-i1) = NC_ObsDiss(See_B03(:,i),See_B05(:,i),...
        Direc_B03(:,i),Direc_B05(:,i),ff,Butm{2},Butm{3},Bdepth{2}(i),...
        Bdepth{3}(i));
    [~,~,~,SeaP.NTED1{2}(:,i+1-i1)] = NC_NeilsonDiss(See_B03(:,i),ff,...
        Bdepth{2}(i),kw,a1,a2,a3,'r');
    [~,~,~,SeaP.NTED2{2}(:,i+1-i1)] = NC_NeilsonDiss(See_B05(:,i),ff,...
        Bdepth{3}(i),kw,a1,a2,a3,'r');
    
    % b/t B05 and B10:
    SeaP.ObsTED{3}(:,i+1-i1) = NC_ObsDiss(See_B05(1:35,i),See_B10(:,i-151),...
        Direc_B05(1:35,i),Direc_B10(:,i),ff(1:35),Butm{3},ADCP.B10.utm,Bdepth{3}(i),...
        Depth_B10(i-151));
    [~,~,~,SeaP.NTED1{3}(:,i+1-i1)] = NC_NeilsonDiss(See_B05(1:35,i),ff(1:35),...
        Bdepth{3}(i),kw,a1,a2,a3,'r');
    [~,~,~,SeaP.NTED2{3}(:,i+1-i1)] = NC_NeilsonDiss(See_B10(:,i-151),ff(1:35),...
        Depth_B10(i-151),kw,a1,a2,a3,'r');
    
    % b/t B10 and B13:
    SeaP.ObsTED{4}(:,i+1-i1) = NC_ObsDiss(See_B10(1:35,i-151),See_B13(:,i-151),...
        Direc_B10(:,i),Direc_B13(:,i),ff(1:35),ADCP.B10.utm,ADCP.B13.utm,Depth_B13(i-151),...
        Depth_B10(i-151));
    [~,~,~,SeaP.NTED1{4}(:,i+1-i1)] = NC_NeilsonDiss(See_B10(:,i-151),ff(1:35),...
        Depth_B10(i-151),kw,a1,a2,a3,'r'); 
    [~,~,~,SeaP.NTED2{4}(:,i+1-i1)] = NC_NeilsonDiss(See_B13(:,i-151),ff(1:35),...
        Depth_B13(i-151),kw,a1,a2,a3,'r'); 
end
SeaP.NTED{1} = (SeaP.NTED1{1} + SeaP.NTED2{1})./2;
SeaP.NTED{2} = (SeaP.NTED1{2} + SeaP.NTED2{2})./2;
SeaP.NTED{3} = (SeaP.NTED1{3} + SeaP.NTED2{3})./2;
SeaP.NTED{4} = (SeaP.NTED1{4} + SeaP.NTED2{4})./2;

    % Find Time Averages:
% From B01 to B03
SeaP.meanObsTED{1} = mean(SeaP.ObsTED{1},2);
SeaP.meanNTED{1} = mean(SeaP.NTED{1},2);
% From B03 to B05
SeaP.meanObsTED{2} = mean(SeaP.ObsTED{2},2);
SeaP.meanNTED{2} = mean(SeaP.NTED{2},2);
% From B05 to B10
SeaP.meanObsTED{3} = mean(SeaP.ObsTED{3},2);
SeaP.meanNTED{3} = mean(SeaP.NTED{3},2);
% From B10 to B13
SeaP.meanObsTED{4} = mean(SeaP.ObsTED{4},2);
SeaP.meanNTED{4} = mean(SeaP.NTED{4},2);

    % Find Stds:
% From B01 to B03
SeaP.stdObsTED{1} = std(SeaP.ObsTED{1},0,2);
SeaP.stdNTED{1} = std(SeaP.NTED{1},0,2);
% From B03 to B05
SeaP.stdObsTED{2} = std(SeaP.ObsTED{2},0,2);
SeaP.stdNTED{2} = std(SeaP.NTED{2},0,2);
% From B05 to B10
SeaP.stdObsTED{3} = std(SeaP.ObsTED{3},0,2);
SeaP.stdNTED{3} = std(SeaP.NTED{3},0,2);
% From B10 to B13
SeaP.stdObsTED{4} = std(SeaP.ObsTED{4},0,2);
SeaP.stdNTED{4} = std(SeaP.NTED{4},0,2);

% Plotting: 
figure(10);clf;
subplot(4,1,1)
plot(ff,SeaP.meanObsTED{1},'r','LineWidth',2)
hold on
plot(ff,SeaP.meanNTED{1},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell')
xlim([0 ff(35)])
ylim([-0.2 0.5])
grid on
title('Sea Peak Period: Dissipation from B01 to B03','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(4,1,2)
plot(ff,SeaP.meanObsTED{2},'r','LineWidth',2)
hold on
plot(ff,SeaP.meanNTED{2},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell')
xlim([0 ff(35)])
ylim([-0.5 1])
grid on
title('Sea Peak Period: Dissipation from B03 to B05','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(4,1,3)
plot(ff(1:35),SeaP.meanObsTED{3},'r','LineWidth',2)
hold on
plot(ff(1:35),SeaP.meanNTED{3},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','southeast')
xlim([0 ff(35)])
ylim([-2 2])
grid on
title('Sea Peak Period: Dissipation from B05 to B10','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(4,1,4)
plot(ff(1:35),SeaP.meanObsTED{4},'r','LineWidth',2)
hold on
plot(ff(1:35),SeaP.meanNTED{4},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','southeast')
xlim([0 ff(35)])
ylim([-2 4])
grid on
title('Sea Peak Period: Dissipation from B10 to B13','fontsize',16)
xlabel('f [Hz]','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)


%% TEDs for Swell Peak Time Period

i1 = SwellP.ind1; % for ADCP this is 91 (- 151)
i2 = SwellP.ind2; % for ADCP this is 130 (- 151)
for i = i1:i2
    % b/t B01 and B03:
    SwellP.ObsTED{1}(:,i+1-i1) = NC_ObsDiss(See_B01_1st(:,i),See_B03(:,i),...
        Direc_B01_1st(:,i),Direc_B03(:,i),ff,Butm{1},Butm{2},Bdepth{1}(i),...
        Bdepth{2}(i));
    [~,~,~,SwellP.NTED1{1}(:,i+1-i1)] = NC_NeilsonDiss(See_B01_1st(:,i),ff,...
        Bdepth{1}(i),kw,a1,a2,a3,'r');
    [~,~,~,SwellP.NTED2{1}(:,i+1-i1)] = NC_NeilsonDiss(See_B03(:,i),ff,...
        Bdepth{2}(i),kw,a1,a2,a3,'r');
    
    % b/t B03 and B05:
    SwellP.ObsTED{2}(:,i+1-i1) = NC_ObsDiss(See_B03(:,i),See_B05(:,i),...
        Direc_B03(:,i),Direc_B05(:,i),ff,Butm{2},Butm{3},Bdepth{2}(i),...
        Bdepth{3}(i));
    [~,~,~,SwellP.NTED1{2}(:,i+1-i1)] = NC_NeilsonDiss(See_B03(:,i),ff,...
        Bdepth{2}(i),kw,a1,a2,a3,'r');
    [~,~,~,SwellP.NTED2{2}(:,i+1-i1)] = NC_NeilsonDiss(See_B05(:,i),ff,...
        Bdepth{3}(i),kw,a1,a2,a3,'r');
    
    % b/t B05 and B10:
    SwellP.ObsTED{3}(:,i+1-i1) = NC_ObsDiss(See_B05(1:35,i),See_B10(:,i-151),...
        Direc_B05(1:35,i),Direc_B10(:,i),ff(1:35),Butm{3},ADCP.B10.utm,Bdepth{3}(i),...
        Depth_B10(i-151));
    [~,~,~,SwellP.NTED1{3}(:,i+1-i1)] = NC_NeilsonDiss(See_B05(1:35,i),ff(1:35),...
        Bdepth{3}(i),kw,a1,a2,a3,'r');
    [~,~,~,SwellP.NTED2{3}(:,i+1-i1)] = NC_NeilsonDiss(See_B10(:,i-151),ff(1:35),...
        Depth_B10(i-151),kw,a1,a2,a3,'r');
    
    % b/t B10 and B13:
    SwellP.ObsTED{4}(:,i+1-i1) = NC_ObsDiss(See_B10(1:35,i-151),See_B13(:,i-151),...
        Direc_B10(:,i),Direc_B13(:,i),ff(1:35),ADCP.B10.utm,ADCP.B13.utm,Depth_B13(i-151),...
        Depth_B10(i-151));
    [~,~,~,SwellP.NTED1{4}(:,i+1-i1)] = NC_NeilsonDiss(See_B10(:,i-151),ff(1:35),...
        Depth_B10(i-151),kw,a1,a2,a3,'r'); 
    [~,~,~,SwellP.NTED2{4}(:,i+1-i1)] = NC_NeilsonDiss(See_B13(:,i-151),ff(1:35),...
        Depth_B13(i-151),kw,a1,a2,a3,'r'); 
end
SwellP.NTED{1} = (SwellP.NTED1{1} + SwellP.NTED2{1})./2;
SwellP.NTED{2} = (SwellP.NTED1{2} + SwellP.NTED2{2})./2;
SwellP.NTED{3} = (SwellP.NTED1{3} + SwellP.NTED2{3})./2;
SwellP.NTED{4} = (SwellP.NTED1{4} + SwellP.NTED2{4})./2;

    % Find Time Averages:
% From B01 to B03
SwellP.meanObsTED{1} = mean(SwellP.ObsTED{1},2);
SwellP.meanNTED{1} = mean(SwellP.NTED{1},2);
% From B03 to B05
SwellP.meanObsTED{2} = mean(SwellP.ObsTED{2},2);
SwellP.meanNTED{2} = mean(SwellP.NTED{2},2);
% From B05 to B10
SwellP.meanObsTED{3} = mean(SwellP.ObsTED{3},2);
SwellP.meanNTED{3} = mean(SwellP.NTED{3},2);
% From B10 to B13
SwellP.meanObsTED{4} = mean(SwellP.ObsTED{4},2);
SwellP.meanNTED{4} = mean(SwellP.NTED{4},2);

    % Find Stds:
% From B01 to B03
SwellP.stdObsTED{1} = std(SwellP.ObsTED{1},0,2);
SwellP.stdNTED{1} = std(SwellP.NTED{1},0,2);
% From B03 to B05
SwellP.stdObsTED{2} = std(SwellP.ObsTED{2},0,2);
SwellP.stdNTED{2} = std(SwellP.NTED{2},0,2);
% From B05 to B10
SwellP.stdObsTED{3} = std(SwellP.ObsTED{3},0,2);
SwellP.stdNTED{3} = std(SwellP.NTED{3},0,2);
% From B10 to B13
SwellP.stdObsTED{4} = std(SwellP.ObsTED{4},0,2);
SwellP.stdNTED{4} = std(SwellP.NTED{4},0,2);

% Plotting: 
figure(11);clf;
subplot(4,1,1)
plot(ff,SwellP.meanObsTED{1},'r','LineWidth',2)
hold on
plot(ff,SwellP.meanNTED{1},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','southeast')
xlim([0 ff(35)])
ylim([-1 0.5])
grid on
title('Swell Peak Period: Dissipation from B01 to B03','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(4,1,2)
plot(ff,SwellP.meanObsTED{2},'r','LineWidth',2)
hold on
plot(ff,SwellP.meanNTED{2},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','southeast')
xlim([0 ff(35)])
ylim([-3.5 1])
grid on
title('Swell Peak Period: Dissipation from B03 to B05','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(4,1,3)
plot(ff(1:35),SwellP.meanObsTED{3},'r','LineWidth',2)
hold on
plot(ff(1:35),SwellP.meanNTED{3},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','southeast')
xlim([0 ff(35)])
ylim([-4 4])
grid on
title('Swell Peak Period: Dissipation from B05 to B10','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(4,1,4)
plot(ff(1:35),SwellP.meanObsTED{4},'r','LineWidth',2)
hold on
plot(ff(1:35),SwellP.meanNTED{4},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','southeast')
xlim([0 ff(35)])
ylim([-5.5 3.5])
grid on
title('Swell Peak Period: Dissipation from B10 to B13','fontsize',16)
xlabel('f [Hz]','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)


%% TEDs for Both Sea and Swell Peak Time Period

i1 = BothP.ind1; % for ADCP this is 207 (- 151)
i2 = BothP.ind2; % for ADCP this is 246 (- 151)
for i = i1:i2
    % b/t B01 and B03:
    BothP.ObsTED{1}(:,i+1-i1) = NC_ObsDiss(See_B01_1st(:,i),See_B03(:,i),...
        Direc_B01_1st(:,i),Direc_B03(:,i),ff,Butm{1},Butm{2},Bdepth{1}(i),...
        Bdepth{2}(i));
    [~,~,~,BothP.NTED1{1}(:,i+1-i1)] = NC_NeilsonDiss(See_B01_1st(:,i),ff,...
        Bdepth{1}(i),kw,a1,a2,a3,'r');
    [~,~,~,BothP.NTED2{1}(:,i+1-i1)] = NC_NeilsonDiss(See_B03(:,i),ff,...
        Bdepth{2}(i),kw,a1,a2,a3,'r');
    
    % b/t B03 and B05:
    BothP.ObsTED{2}(:,i+1-i1) = NC_ObsDiss(See_B03(:,i),See_B05(:,i),...
        Direc_B03(:,i),Direc_B05(:,i),ff,Butm{2},Butm{3},Bdepth{2}(i),...
        Bdepth{3}(i));
    [~,~,~,BothP.NTED1{2}(:,i+1-i1)] = NC_NeilsonDiss(See_B03(:,i),ff,...
        Bdepth{2}(i),kw,a1,a2,a3,'r');
    [~,~,~,BothP.NTED2{2}(:,i+1-i1)] = NC_NeilsonDiss(See_B05(:,i),ff,...
        Bdepth{3}(i),kw,a1,a2,a3,'r');
    
    % b/t B05 and B10:
    BothP.ObsTED{3}(:,i+1-i1) = NC_ObsDiss(See_B05(1:35,i),See_B10(:,i-151),...
        Direc_B05(1:35,i),Direc_B10(:,i),ff(1:35),Butm{3},ADCP.B10.utm,Bdepth{3}(i),...
        Depth_B10(i-151));
    [~,~,~,BothP.NTED1{3}(:,i+1-i1)] = NC_NeilsonDiss(See_B05(1:35,i),ff(1:35),...
        Bdepth{3}(i),kw,a1,a2,a3,'r');
    [~,~,~,BothP.NTED2{3}(:,i+1-i1)] = NC_NeilsonDiss(See_B10(:,i-151),ff(1:35),...
        Depth_B10(i-151),kw,a1,a2,a3,'r');
    
    % b/t B10 and B13:
    BothP.ObsTED{4}(:,i+1-i1) = NC_ObsDiss(See_B10(1:35,i-151),See_B13(:,i-151),...
        Direc_B10(:,i),Direc_B13(:,i),ff(1:35),ADCP.B10.utm,ADCP.B13.utm,Depth_B13(i-151),...
        Depth_B10(i-151));
    [~,~,~,BothP.NTED1{4}(:,i+1-i1)] = NC_NeilsonDiss(See_B10(:,i-151),ff(1:35),...
        Depth_B10(i-151),kw,a1,a2,a3,'r'); 
    [~,~,~,BothP.NTED2{4}(:,i+1-i1)] = NC_NeilsonDiss(See_B13(:,i-151),ff(1:35),...
        Depth_B13(i-151),kw,a1,a2,a3,'r'); 
end
BothP.NTED{1} = (BothP.NTED1{1} + BothP.NTED2{1})./2;
BothP.NTED{2} = (BothP.NTED1{2} + BothP.NTED2{2})./2;
BothP.NTED{3} = (BothP.NTED1{3} + BothP.NTED2{3})./2;
BothP.NTED{4} = (BothP.NTED1{4} + BothP.NTED2{4})./2;

    % Find Time Averages:
% From B01 to B03
BothP.meanObsTED{1} = mean(BothP.ObsTED{1},2);
BothP.meanNTED{1} = mean(BothP.NTED{1},2);
% From B03 to B05
BothP.meanObsTED{2} = mean(BothP.ObsTED{2},2);
BothP.meanNTED{2} = mean(BothP.NTED{2},2);
% From B05 to B10
BothP.meanObsTED{3} = mean(BothP.ObsTED{3},2);
BothP.meanNTED{3} = mean(BothP.NTED{3},2);
% From B10 to B13
BothP.meanObsTED{4} = mean(BothP.ObsTED{4},2);
BothP.meanNTED{4} = mean(BothP.NTED{4},2);

    % Find Stds:
% From B01 to B03
BothP.stdObsTED{1} = std(BothP.ObsTED{1},0,2);
BothP.stdNTED{1} = std(BothP.NTED{1},0,2);
% From B03 to B05
BothP.stdObsTED{2} = std(BothP.ObsTED{2},0,2);
BothP.stdNTED{2} = std(BothP.NTED{2},0,2);
% From B05 to B10
BothP.stdObsTED{3} = std(BothP.ObsTED{3},0,2);
BothP.stdNTED{3} = std(BothP.NTED{3},0,2);
% From B10 to B13
BothP.stdObsTED{4} = std(BothP.ObsTED{4},0,2);
BothP.stdNTED{4} = std(BothP.NTED{4},0,2);

% Plotting: 
figure(12);clf;
subplot(4,1,1)
plot(ff,BothP.meanObsTED{1},'r','LineWidth',2)
hold on
plot(ff,BothP.meanNTED{1},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell')
xlim([0 ff(35)])
ylim([-1 0.75])
grid on
title('Sea/Swell Peak Period: Dissipation from B01 to B03','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(4,1,2)
plot(ff,BothP.meanObsTED{2},'r','LineWidth',2)
hold on
plot(ff,BothP.meanNTED{2},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','northeast')
xlim([0 ff(35)])
ylim([-1.5 2])
grid on
title('Sea/Swell Peak Period: Dissipation from B03 to B05','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(4,1,3)
plot(ff(1:35),BothP.meanObsTED{3},'r','LineWidth',2)
hold on
plot(ff(1:35),BothP.meanNTED{3},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','northeast')
xlim([0 ff(35)])
ylim([-6 3.5])
grid on
title('Sea/Swell Peak Period: Dissipation from B05 to B10','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(4,1,4)
plot(ff(1:35),BothP.meanObsTED{4},'r','LineWidth',2)
hold on
plot(ff(1:35),BothP.meanNTED{4},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','northeast')
xlim([0 ff(35)])
ylim([-4 14])
grid on
title('Sea/Swell Peak Period: Dissipation from B10 to B13','fontsize',16)
xlabel('f [Hz]','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)


%% Save Variables and Figures

clear A an an1 an2 Bothx1 Bothx2 Seax1 Seax2 Swellx1 Swellx2 spotterL2...
    cb i i1 i2 SM yL

%     % Save Variables
% save('3Periods.mat','SeaP','SwellP','BothP')
% 
%     % Save Figures
% cd 'Spring Figures/ErrBarAvgSpec_PLOTS'
% % Spectrogram Figure
% saveas(figure(20),'Periods_on_Spectro.jpeg')
% % Sea Period Figures
% saveas(figure(1),'B01_SeaP.jpeg')
% saveas(figure(4),'B03_SeaP.jpeg')
% saveas(figure(7),'B05_SeaP.jpeg')
% % Swell Period Figures
% saveas(figure(2),'B01_SwellP.jpeg')
% saveas(figure(5),'B03_SwellP.jpeg')
% saveas(figure(8),'B05_SwellP.jpeg')
% % Both Period Figures
% saveas(figure(3),'B01_BothP.jpeg')
% saveas(figure(6),'B03_BothP.jpeg')
% saveas(figure(9),'B05_BothP.jpeg')
% % Dissipation Figures
% saveas(figure(10),'SeaP_Diss.jpeg')
% saveas(figure(11),'SwellP_Diss.jpeg')
% saveas(figure(12),'BothP_Diss.jpeg')









