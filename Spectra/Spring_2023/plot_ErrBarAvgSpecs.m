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
% TO-DO: 
%   - add err bars to the TED plots
%   - investigate these TED plots
%
% Questions/Notes:
%           - How should I add err bars to the TED plots without them being
%             too crowded?

clc;clear;

%% Suanda Example (Basic)
% original example of envelope function (bgpatch2)

addpath('../Summer (0516-0728)/Start Data')
load('roxsi_spotter_L2_B05_1153_reduced.mat')

% Calculate the time-average and standard error of the energy:
N = length(spotterL2.dtime);
mn = mean(spotterL2.See,2);
st = std(spotterL2.See,0,2)/N;

% Plot the mean and standard error envelope:
figure(1);clf;
plot(spotterL2.frequency,mn,'k')
hold on 
bgpatch2(spotterL2.frequency,mn,st,'r')

close 
clear mn N st


%% Preliminaries
% The preliminaries for the original plots of the 3 selected time periods

% Load Normal Wave Direction: 
load('WBvariables.mat','BNormWaveDir')

    % Define the difference between Sea and Swell:
Def_Freq = 0.08; % >: Sea       % <: Swell

    % Fully Time Averaged Variables (For B05):
FTA_See = mean(spotterL2.See,2);
FTA_Dir = meanangle(spotterL2.EMEM.meandir_f,2) + 360;
FTA_DirSpread = meanangle(spotterL2.EMEM.meandirspread_f,2);

    % Frequency Vector: 
ff = spotterL2.frequency;


%% B05: Sea Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy only peaks in the sea band
%
% - Time: June 17 @ 18:00 - June 19 @ 09:00   (1 day and 15 hours)


    % Time Indexes:
SeaP.ind1 = 55;
SeaP.ind2 = 94;
SeaP.N = 94 - 55 + 1;
SeaP.TimeStart = spotterL2.dtime(SeaP.ind1); % June 17 @ 18:00
SeaP.TimeEnd = spotterL2.dtime(SeaP.ind2); % June 19 @ 09:00

    % Assign Variables (for the day and a half):
% Energy
SeaP.mSee = mean(spotterL2.See(:,SeaP.ind1:SeaP.ind2),2);
SeaP.std_mSee = std(spotterL2.See(:,SeaP.ind1:SeaP.ind2),0,2);
% Direction
SeaP.mDir = meanangle(spotterL2.EMEM.meandir_f(:,SeaP.ind1:SeaP.ind2),2) + 360;
SeaP.std_mDir = std(spotterL2.EMEM.meandir_f(:,SeaP.ind1:SeaP.ind2),0,2);
% Directional Spread
SeaP.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,SeaP.ind1:SeaP.ind2),2);
SeaP.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,SeaP.ind1:SeaP.ind2),0,2);

    % Plot:
figure(2);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,SeaP.mSee,SeaP.std_mSee,'k2')
hold on
plot(ff,SeaP.mSee,'r','LineWidth',1.5)
plot(ff,FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
grid on
title('B05 Sea Peak Period: Mean Energy (Jun 17 @ 18:00 - Jun 19 @ 09:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.5 4])
xlim([0 ff(37)])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,SeaP.mDir,SeaP.std_mDir,'k2')
hold on
plot(ff,SeaP.mDir,'b','LineWidth',1.5)
plot(ff,FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
grid on
title('B05 Sea Peak Period: Mean Direction (Jun 17 @ 18:00 - Jun 19 @ 09:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([190 325])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','south')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,SeaP.mDirSpread,SeaP.std_mDirSpread,'k2')
hold on
plot(ff,SeaP.mDirSpread,'g','LineWidth',1.5)
plot(ff,FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
grid on
title('B05 Sea Peak Period: Mean Directional Spread (Jun 17 @ 18:00 - Jun 19 @ 09:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (17th-19th)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% B05: Swell Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy only peaks in the swell band
%
% - Time: June 25 @ 13:00 - June 27 @ 04:00      (1 day and 15 hours)

    % Time Indexes:
SwellP.ind1 = 242;
SwellP.ind2 = 281;
SwellP.N = 94 - 55 + 1;
SwellP.TimeStart = spotterL2.dtime(SwellP.ind1); % June 25 @ 13:00
SwellP.TimeEnd = spotterL2.dtime(SwellP.ind2); % June 27 @ 04:00

    % Assign Variables (for the day and a half):
% Energy
SwellP.mSee = mean(spotterL2.See(:,SwellP.ind1:SwellP.ind2),2);
SwellP.std_mSee = std(spotterL2.See(:,SwellP.ind1:SwellP.ind2),0,2);
% Direction
SwellP.mDir = meanangle(spotterL2.EMEM.meandir_f(:,SwellP.ind1:SwellP.ind2),2) + 360;
SwellP.std_mDir = std(spotterL2.EMEM.meandir_f(:,SwellP.ind1:SwellP.ind2),0,2);
% Directional Spread
SwellP.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,SwellP.ind1:SwellP.ind2),2);
SwellP.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,SwellP.ind1:SwellP.ind2),0,2);

    % Plot:
figure(3);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,SwellP.mSee,SwellP.std_mSee,'k2')
hold on
plot(ff,SwellP.mSee,'r','LineWidth',1.5)
plot(ff,FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
grid on
title('B05 Swell Peak Period: Mean Energy (June 25 @ 13:00 - June 27 @ 04:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.25 1.25])
xlim([0 ff(37)])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,SwellP.mDir,SwellP.std_mDir,'k2')
hold on
plot(ff,SwellP.mDir,'b','LineWidth',1.5)
plot(ff,FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
grid on
title('B05 Swell Peak Period: Mean Direction (June 25 @ 13:00 - June 27 @ 04:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([170 350])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','south')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,SwellP.mDirSpread,SwellP.std_mDirSpread,'k2')
hold on
plot(ff,SwellP.mDirSpread,'g','LineWidth',1.5)
plot(ff,FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
grid on
title('B05 Swell Peak Period: Mean Directional Spread (June 25 @ 13:00 - June 27 @ 04:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (25th-27th)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% B05: Swell and Sea Band Peak Time
% - plot the Energy, Direction, and Directional Spread for the time where
%   energy peaks in both the swell and the sea bands
%
% - Time: Jun 30 @ 09:00 - July 02 @ 00:00      (1 day and 15 hours)

    % Time Indexes:
BothP.ind1 = 358;
BothP.ind2 = 397;
BothP.N = 94 - 55 + 1;
BothP.TimeStart = spotterL2.dtime(BothP.ind1); % June 30 @ 09:00
BothP.TimeEnd = spotterL2.dtime(BothP.ind2); % July 02 @ 00:00

    % Assign Variables (for the day and a half):
% Energy
BothP.mSee = mean(spotterL2.See(:,BothP.ind1:BothP.ind2),2);
BothP.std_mSee = std(spotterL2.See(:,BothP.ind1:BothP.ind2),0,2);
% Direction
BothP.mDir = meanangle(spotterL2.EMEM.meandir_f(:,BothP.ind1:BothP.ind2),2) + 360;
BothP.std_mDir = std(spotterL2.EMEM.meandir_f(:,BothP.ind1:BothP.ind2),0,2);
% Directional Spread
BothP.mDirSpread = meanangle(spotterL2.EMEM.meandirspread_f(:,BothP.ind1:BothP.ind2),2);
BothP.std_mDirSpread = std(spotterL2.EMEM.meandirspread_f(:,BothP.ind1:BothP.ind2),0,2);

    % Plot:
figure(4);clf;
%Energy
subplot(3,1,1)
bgpatch2(ff,BothP.mSee,BothP.std_mSee,'k2')
hold on
plot(ff,BothP.mSee,'r','LineWidth',1.5)
plot(ff,FTA_See,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
grid on
title('B05 Both Peaks Period: Mean Energy (Jun 30 @ 09:00 - Jul 02 @ 00:00)','fontsize',16)
ylabel('E (m^2/Hz)','fontsize',16)
ylim([-0.25 1.6])
xlim([0 ff(37)])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell')
% Direction
subplot(3,1,2)
bgpatch2(ff,BothP.mDir,BothP.std_mDir,'k2')
hold on
plot(ff,BothP.mDir,'b','LineWidth',1.5)
plot(ff,FTA_Dir,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
yline(BNormWaveDir,'-k')
grid on
title('B05 Both Peaks Period: Mean Direction (Jun 30 @ 09:00 - Jul 02 @ 00:00)','fontsize',16)
ylabel('Dir (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([180 360])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell','Shore Normal',...
    'location','south')
% Directional Spread
subplot(3,1,3)
bgpatch2(ff,BothP.mDirSpread,BothP.std_mDirSpread,'k2')
hold on
plot(ff,BothP.mDirSpread,'g','LineWidth',1.5)
plot(ff,FTA_DirSpread,'--k','LineWidth',1.5)
xline(Def_Freq,'--m')
grid on
title('B05 Both Peaks Period: Mean Directional Spread (Jun 30 @ 09:00 - Jul 02 @ 00:00)','fontsize',16)
xlabel('f (Hz)','fontsize',16)
ylabel('Dir Spread (^o)','fontsize',16)
xlim([0 ff(37)])
ylim([18 85])
legend('Avg (30th-2nd)','Avg (34 Days)','Sea vs Swell','location','southeast')


%% Preliminaries (for TED and Nielsen TED)

% Add the path for the functions
addpath('../Summer (0516-0728)/Main Functions')

% Load data (1st set) for B01
load('roxsi_spotter_L2_B01_1158_reduced.mat')
See_B01_1st = spotterL2.See;
Direc_B01_1st = spotterL2.EMEM.meandir_f;
% Load data (2nd set) for B01
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

% Load data needed for functions
load('WBvariables.mat','Butm','Bdepth')

% More needed for functions
kw = 1;   a1 = 5.5;   a2 = -0.2;    a3 = -6.3;

% Have to change ff from starting at 0 to starting at 0.0098 or else the
%  Nielsen function will return all NaNs:
ff = 0.0098:0.0098:(0.0098*129);


%% TEDs for Sea Peak Time Period

i1 = SeaP.ind1;
i2 = SeaP.ind2;
for i = i1:i2
    % b/t B01 and B03:
    SeaP.ObsTED{1}(:,i+1-i1) = NC_ObsDiss(See_B01_1st(:,i),See_B03(:,i),...
        Direc_B01_1st(:,i),Direc_B03(:,i),ff,Butm{1},Butm{2},Bdepth{1}(i),...
        Bdepth{2}(i));
    [~,~,~,SeaP.NTED1{1}(:,i+1-i1)] = NC_NeilsonDiss(See_B01_1st(:,i),ff,...
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
end
SeaP.NTED{1} = (SeaP.NTED1{1} + SeaP.NTED2{1})./2;
SeaP.NTED{2} = (SeaP.NTED1{2} + SeaP.NTED2{2})./2;

    % Find Time Averages:
% From B01 to B03
SeaP.meanObsTED{1} = mean(SeaP.ObsTED{1},2);
SeaP.meanNTED{1} = mean(SeaP.NTED{1},2);
% From B03 to B05
SeaP.meanObsTED{2} = mean(SeaP.ObsTED{2},2);
SeaP.meanNTED{2} = mean(SeaP.NTED{2},2);

    % Find Stds:
% From B01 to B03
SeaP.stdObsTED{1} = std(SeaP.ObsTED{1},0,2);
SeaP.stdNTED{1} = std(SeaP.NTED{1},0,2);
% From B03 to B05
SeaP.stdObsTED{2} = std(SeaP.ObsTED{2},0,2);
SeaP.stdNTED{2} = std(SeaP.NTED{2},0,2);

% Plotting: 
figure(5);clf;
subplot(2,1,1)
plot(ff,SeaP.meanObsTED{1},'r','LineWidth',2)
hold on
plot(ff,SeaP.meanNTED{1},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell')
xlim([0 ff(35)])
ylim([-0.30 1.5])
grid on
title('Sea Peak Period: Dissipation from B01 to B03','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(2,1,2)
plot(ff,SeaP.meanObsTED{2},'r','LineWidth',2)
hold on
plot(ff,SeaP.meanNTED{2},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell')
xlim([0 ff(35)])
ylim([-1 3.5])
grid on
title('Sea Peak Period: Dissipation from B03 to B05','fontsize',16)
xlabel('f [Hz]','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)



%% TEDs for Swell Peak Time Period

i1 = SwellP.ind1;
i2 = SwellP.ind2;
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
end
SwellP.NTED{1} = (SwellP.NTED1{1} + SwellP.NTED2{1})./2;
SwellP.NTED{2} = (SwellP.NTED1{2} + SwellP.NTED2{2})./2;

    % Find Time Averages:
% From B01 to B03
SwellP.meanObsTED{1} = mean(SwellP.ObsTED{1},2);
SwellP.meanNTED{1} = mean(SwellP.NTED{1},2);
% From B03 to B05
SwellP.meanObsTED{2} = mean(SwellP.ObsTED{2},2);
SwellP.meanNTED{2} = mean(SwellP.NTED{2},2);

    % Find Stds:
% From B01 to B03
SwellP.stdObsTED{1} = std(SwellP.ObsTED{1},0,2);
SwellP.stdNTED{1} = std(SwellP.NTED{1},0,2);
% From B03 to B05
SwellP.stdObsTED{2} = std(SwellP.ObsTED{2},0,2);
SwellP.stdNTED{2} = std(SwellP.NTED{2},0,2);

% Plotting: 
figure(6);clf;
subplot(2,1,1)
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

subplot(2,1,2)
plot(ff,SwellP.meanObsTED{2},'r','LineWidth',2)
hold on
plot(ff,SwellP.meanNTED{2},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','southeast')
xlim([0 ff(35)])
ylim([-3.5 1])
grid on
title('Swell Peak Period: Dissipation from B03 to B05','fontsize',16)
xlabel('f [Hz]','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)


%% TEDs for Both Sea and Swell Peak Time Period

i1 = BothP.ind1;
i2 = BothP.ind2;
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
end
BothP.NTED{1} = (BothP.NTED1{1} + BothP.NTED2{1})./2;
BothP.NTED{2} = (BothP.NTED1{2} + BothP.NTED2{2})./2;

    % Find Time Averages:
% From B01 to B03
BothP.meanObsTED{1} = mean(BothP.ObsTED{1},2);
BothP.meanNTED{1} = mean(BothP.NTED{1},2);
% From B03 to B05
BothP.meanObsTED{2} = mean(BothP.ObsTED{2},2);
BothP.meanNTED{2} = mean(BothP.NTED{2},2);

    % Find Stds:
% From B01 to B03
BothP.stdObsTED{1} = std(BothP.ObsTED{1},0,2);
BothP.stdNTED{1} = std(BothP.NTED{1},0,2);
% From B03 to B05
BothP.stdObsTED{2} = std(BothP.ObsTED{2},0,2);
BothP.stdNTED{2} = std(BothP.NTED{2},0,2);

% Plotting: 
figure(7);clf;
subplot(2,1,1)
plot(ff,BothP.meanObsTED{1},'r','LineWidth',2)
hold on
plot(ff,BothP.meanNTED{1},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell')
xlim([0 ff(35)])
ylim([-1 1])
grid on
title('Sea/Swell Peak Period: Dissipation from B01 to B03','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)

subplot(2,1,2)
plot(ff,BothP.meanObsTED{2},'r','LineWidth',2)
hold on
plot(ff,BothP.meanNTED{2},'b','LineWidth',2)
xline(Def_Freq,'--m')
legend('Observed', 'Nielsen','Sea vs Swell','location','northeast')
xlim([0 ff(35)])
ylim([-1.25 1.5])
grid on
title('Sea/Swell Peak Period: Dissipation from B03 to B05','fontsize',16)
xlabel('f [Hz]','fontsize',16)
ylabel('E [m^2/Hz]','fontsize',16)
















