%% DirSpread.m
%
% Noah Clark        12/15/23
%
%
% Purpose: Determine and plot the time-averaged mean directional spread 
%          for all spotter buoys
%
% Note: I don't think that I need to use the meanangle function for these
%       because none of the angles seem to bounce between 0 and 360
%

%% Preliminaries

clc;clear;

addpath('/Users/noahclark/Library/CloudStorage/OneDrive-UNC-Wilmington/CMS Work/Summer (0516-0728)/Start Data')

%% China Rock Buoys

    % B01
% First half of time:
load('roxsi_spotter_L2_B01_1150_reduced.mat')
ff = (1:129).*spotterL2.df; % general for all buoys
B01Spread1 = spotterL2.EMEM.meandirspread_f;
% Second half of time:
load('roxsi_spotter_L2_B01_1158_reduced.mat')
B01Spread2 = spotterL2.EMEM.meandirspread_f;
% Concatonate the two halfs:
B01Spread = cat(2,B01Spread1,B01Spread2);
% Average: 
MeanSpread.B01 = mean(B01Spread,2);

    % B03
load('roxsi_spotter_L2_B03_1152_reduced.mat')
MeanSpread.B03 = mean(spotterL2.EMEM.meandirspread_f,2);

    % B05
load('roxsi_spotter_L2_B05_1153_reduced.mat')
MeanSpread.B05 = mean(spotterL2.EMEM.meandirspread_f,2);


    % Plot the 3 China Rocky buoy mean spreads:
figure(1);clf;
plot(ff,MeanSpread.B01,'r','LineWidth',2)
hold on; grid on;
plot(ff,MeanSpread.B03,'b','LineWidth',2)
plot(ff,MeanSpread.B05,'g','LineWidth',2)
xlabel('Frequency (Hz)','fontsize',17)
ylabel('Dir Spread (^o)','fontsize',15)
title('China Rock Mean Directional Spread','fontsize',20)
legend('B01','B03','B05','fontsize',15)
xlim([0 0.35])


%% Asilomar Buoys

    % X01
load('roxsi_spotter_L2_X01_1151_reduced.mat')
MeanSpread.X01 = mean(spotterL2.EMEM.meandirspread_f,2);

    % X03
load('roxsi_spotter_L2_X03_1157_reduced.mat')
MeanSpread.X03 = nanmean(spotterL2.EMEM.meandirspread_f,2);

    % X04
load('roxsi_spotter_L2_X04_1155_reduced.mat')
MeanSpread.X04 = mean(spotterL2.EMEM.meandirspread_f,2);


    % Plot the 3 Asilomar buoy mean spreads:
figure(2);clf;
plot(ff,MeanSpread.X01,'r','LineWidth',2)
hold on; grid on;
plot(ff,MeanSpread.X03,'b','LineWidth',2)
plot(ff,MeanSpread.X04,'g','LineWidth',2)
xlabel('Frequency (Hz)','fontsize',17)
ylabel('Dir Spread (^o)','fontsize',15)
title('Asilomar Mean Directional Spread','fontsize',20)
legend('X01','X03','X04','fontsize',15)
xlim([0 0.35])


clear B01Spread B01Spread1 B01Spread2 spotterL2


%% Save Variables

save('WBvariables.mat','-append','MeanSpread')



