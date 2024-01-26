%% plot_meanDir.m
%
% Noah Clark
% Date Created: 9/22/23
%
% Purpose: To determine and plot the mean wave direction spectrums at
%           buoy B01 and X01


%% Preliminaries

clc;clear;
ff = 0.0098:0.0098:1.2642;
load('WBvariables.mat','BEMEM','XEMEM')


%% Calculations

for i = 1:3
    BEMEMmean{i} = meanangle(BEMEM{i}');
    if i == 2
        XEMEMmean{i} = meanangle(XEMEM{i}(:,[1:172 174:end])');
    else
        XEMEMmean{i} = meanangle(XEMEM{i}');
    end
    for j = 1:129
        if BEMEMmean{i}(j) < 0
            BEMEMmean{i}(j) = BEMEMmean{i}(j) + 360;
        end
        if XEMEMmean{i}(j) < 0 
            XEMEMmean{i}(j) = XEMEMmean{i}(j) + 360;
        end
    end
end


%% Plotting

figure(1);clf;
plot(ff,BEMEMmean{1},'Linewidth',2)
title('Time Averaged Spectrum for B01')
xlabel('Frequency (Hz)')
ylabel('Direction (degrees)')
grid on
ylim([50 370])
xlim([0 1.3])

figure(2);clf;
plot(ff,XEMEMmean{1},'r','Linewidth',2)
title('Time Averaged Spectrum for X01')
xlabel('Frequency (Hz)')
ylabel('Direction (degrees)')
grid on
ylim([50 370])
xlim([0 1.3])



save('WBvariables.mat','-append','BEMEMmean','XEMEMmean')