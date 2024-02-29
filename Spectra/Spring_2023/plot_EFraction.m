%% plot_EFraction.m
%
% Noah Clark            2/2/2024
%
%
% Purpose: To calculate and plot the fraction of total energy that is found
%          at each day in the sea and swell bands. This is done only for
%          buoy B03.
%
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

%% Preliminaries

clc;clear;

load('WBvariables.mat','BSee','Btime')

See = BSee{2};
Time = Btime{2};

clear Btime

ff = 0:0.0098:(0.0098*128);

%% Calculate

% - The frequency indecies for Sea: 0.088 - 0.245 Hz       (ff(10)-ff(25))
% - The frequency indecies for Swell: 0.0392 - 0.088 Hz     (ff(5)-ff(10))

SeaSee = See(10:25,:);
SwellSee = See(5:10,:);

% Sum all frequencies per hour
fSumSee = sum(See,1);
fSumSeaSee = sum(SeaSee,1);
fSumSwellSee = sum(SwellSee,1);

% Sum Daily
for i = 1:34
    dSumSee(i) = sum(fSumSee(i*24-23:i*24));
    dSumSeaSee(i) = sum(fSumSeaSee(i*24-23:i*24));
    dSumSwellSee(i) = sum(fSumSwellSee(i*24-23:i*24));
end

% Calulate Fraction
FracSeaSee = dSumSeaSee./dSumSee;
FracSwellSee = dSumSwellSee./dSumSee;
% Calculate the Average Fractions (from these daily averages)
avgFracSeaSee = mean(FracSeaSee);
avgFracSwellSee = mean(FracSwellSee);

% temporary time vector
Tday = datetime(2022,6,16):datetime(2022,7,19);


%% Plot

figure(1);clf;
plot(Tday,FracSeaSee,'xk','MarkerSize',10,'LineWidth',2)
hold on; grid on
plot(Tday,FracSwellSee,'xr','MarkerSize',10,'LineWidth',2)
legend('Sea','Swell','fontsize',14,'location','east')
title('B03: Energy Fraction in Sea & Swell','fontsize',20,'fontweight','normal')
xlabel('Date','fontsize',18)
ylabel('Energy Fraction (E_s_e_a_/_s_w_e_l_l/E_t_o_t)','fontsize',16)
xlim([datetime(2022,6,14) datetime(2022,7,21)])


%% Calculate the average percentages (not daily avgs)

% B01
B01FracSea = sum(BSee{1}(10:25,:),'all')/sum(BSee{1},'all'); % 0.792
B01FracSwell = sum(BSee{1}(5:10,:),'all')/sum(BSee{1},'all'); % 0.104
% B03
B03FracSea = sum(BSee{2}(10:25,:),'all')/sum(BSee{2},'all'); % 0.759
B03FracSwell = sum(BSee{2}(5:10,:),'all')/sum(BSee{2},'all'); % 0.130
% B05
B05FracSea = sum(BSee{3}(10:25,:),'all')/sum(BSee{2},'all'); % 0.703
B05FracSwell = sum(BSee{3}(5:10,:),'all')/sum(BSee{2},'all'); % 0.142

% B10
addpath('../Fall (0825)/ADCP Data')
load('roxsi_signature_L2_B10_See.mat')
See_B10 = A.See;
B10FracSea = nansum(See_B10(33:89,:),'all')/nansum(See_B10,'all') % 0.748
B10FracSwell = nansum(See_B10(15:33,:),'all')/nansum(See_B10,'all') % 0.112
% B13
load('roxsi_signature_L2_B13_See.mat')
See_B13 = A.See;
B13FracSea = nansum(See_B13(33:89,:),'all')/nansum(See_B13,'all') % 0.687
B13FracSwell = nansum(See_B13(15:33,:),'all')/nansum(See_B13,'all') % 0.206

% Average for all 5
AVGFracSea = mean([B01FracSea B03FracSea B05FracSea B10FracSea B13FracSea]); % 0.738
AVGFracSwell = mean([B01FracSwell B03FracSwell B05FracSwell B10FracSwell B13FracSwell]) % 0.139

%% Save Figure

cd 'Spring Figures/Daily_Averages'
saveas(figure(1),'TotalE_fraction.jpeg')








