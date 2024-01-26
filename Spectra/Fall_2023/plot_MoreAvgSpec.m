%% plot_MoreAvgSpec.m
%
% Noah Clark
%
% Purpose: - Create time-averaged spectrums for all spotter buoys and ADCPs
%          - In one figure (1 for CR and 1 for A) plot time-averaged
%             directions 


clc;clear;


% To-Do: add a line for normal to shoreline for the wind plots

%% Data Organization

%addpath('C:\Users\nsc4750\OneDrive - UNC-Wilmington\CMS Work\Summer (0516-0728)')
load('WBvariables.mat','Bmeanspec','Xmeanspec','Bfreq','Xfreq',...
    'BEMEMmean','XEMEMmean','BNormWaveDir','XNormWaveDir')

addpath('../Fall (0825)/ADCP Data')
ADCPnames = {'roxsi_signature_L2_B10_See.mat',...
    'roxsi_signature_L2_B13_See.mat','roxsi_signature_L2_X05_See.mat',...
    'roxsi_signature_L2_X11_See.mat'};
ADCPnamesR = {'B10', 'B13', 'X05', 'X11'};


addpath('../Fall (0825)/SmartMooring Data')
SMnames = {'roxsi_smartmooring_L2_E01sp_1851.mat','roxsi_smartmooring_L2_E02sp_1859.mat',...
    'roxsi_smartmooring_L2_E05sp_1853.mat','roxsi_smartmooring_L2_E07sp_1855.mat',...
    'roxsi_smartmooring_L2_E07sp_1857.mat','roxsi_smartmooring_L2_E08sp_1852.mat',...
    'roxsi_smartmooring_L2_E09sp_1850.mat','roxsi_smartmooring_L2_E09sp_1856.mat',...
    'roxsi_smartmooring_L2_E10sp_1848.mat','roxsi_smartmooring_L2_E11sp_1860.mat',...
    'roxsi_smartmooring_L2_E13sp_1849.mat'};
SMnamesR = {'E01','E02','E05','E07','E08','E09','E10','E11','E13'};

proj = projcrs(26944);

    % Smart Moorings:
load(SMnames{1});
SM.E01.See = spotsmartL2.Spp;
SMt = spotsmartL2.dtime;
SMf = spotsmartL2.frequency;
SM.time = SMt;
SM.freq = SMf;
SM.E01.Lat = spotsmartL2.latitude;
SM.E01.Lon = spotsmartL2.longitude;
[xutm,yutm] = projfwd(proj,SM.E01.Lat,SM.E01.Lon);
SM.E01.utm = [xutm,yutm];
SM.E01.depth = spotsmartL2.mean_depth;

load(SMnames{2});
SM.E02.See = spotsmartL2.Spp;
SM.E02.Lat = spotsmartL2.latitude;
SM.E02.Lon = spotsmartL2.longitude;
[xutm,yutm] = projfwd(proj,SM.E02.Lat,SM.E02.Lon);
SM.E02.utm = [xutm,yutm];
SM.E02.depth = spotsmartL2.mean_depth;

load(SMnames{3});
SM.E05.See = spotsmartL2.Spp;
SM.E05.Lat = spotsmartL2.latitude;
SM.E05.Lon = spotsmartL2.longitude;
[xutm,yutm] = projfwd(proj,SM.E05.Lat,SM.E05.Lon);
SM.E05.utm = [xutm,yutm];
SM.E05.depth = spotsmartL2.mean_depth;

load(SMnames{4});
SMtemp1 = spotsmartL2;
load(SMnames{5});
SMtemp2 = spotsmartL2;
SM.E07.See = cat(2,SMtemp1.Spp(:,1:202),SMtemp2.Spp(:,203:end));
SM.E07.Lat = spotsmartL2.latitude;
SM.E07.Lon = spotsmartL2.longitude;
[xutm,yutm] = projfwd(proj,SM.E07.Lat,SM.E07.Lon);
SM.E07.utm = [xutm,yutm];
SM.E07.depth = spotsmartL2.mean_depth;

load(SMnames{6});
SM.E08.See = spotsmartL2.Spp;
SM.E08.Lat = spotsmartL2.latitude;
SM.E08.Lon = spotsmartL2.longitude;
[xutm,yutm] = projfwd(proj,SM.E08.Lat,SM.E08.Lon);
SM.E08.utm = [xutm,yutm];
SM.E08.depth = spotsmartL2.mean_depth;

load(SMnames{7});
SMtemp1 = spotsmartL2;
load(SMnames{8});
SMtemp2 = spotsmartL2;
SM.E09.See = cat(2,SMtemp1.Spp(:,1:92),SMtemp2.Spp(:,93:end));
SM.E09.Lat = spotsmartL2.latitude;
SM.E09.Lon = spotsmartL2.longitude;
[xutm,yutm] = projfwd(proj,SM.E09.Lat,SM.E09.Lon);
SM.E09.utm = [xutm,yutm];
SM.E09.depth = spotsmartL2.mean_depth;

load(SMnames{9});
SM.E10.See = spotsmartL2.Spp;
SM.E10.Lat = spotsmartL2.latitude;
SM.E10.Lon = spotsmartL2.longitude;
[xutm,yutm] = projfwd(proj,SM.E10.Lat,SM.E10.Lon);
SM.E10.utm = [xutm,yutm];
SM.E10.depth = spotsmartL2.mean_depth;

load(SMnames{10});
SM.E11.See = spotsmartL2.Spp;
SM.E11.Lat = spotsmartL2.latitude;
SM.E11.Lon = spotsmartL2.longitude;
[xutm,yutm] = projfwd(proj,SM.E11.Lat,SM.E11.Lon);
SM.E11.utm = [xutm,yutm];
SM.E11.depth = spotsmartL2.mean_depth;

load(SMnames{11});
SM.E13.See = spotsmartL2.Spp;
SM.E13.Lat = spotsmartL2.latitude;
SM.E13.Lon = spotsmartL2.longitude;
[xutm,yutm] = projfwd(proj,SM.E13.Lat,SM.E13.Lon);
SM.E13.utm = [xutm,yutm];
SM.E13.depth = spotsmartL2.mean_depth;


    % ADCPs:
load(ADCPnames{1});
ADCPt = A.time;
ADCPf = A.freq;  
ADCP.freq = ADCPf;
ADCP.time = ADCPt;
ADCP.B10.See = A.See;
ADCP.B10.Lat = A.latitude;
ADCP.B10.Lon = A.longitude;
[xutm,yutm] = projfwd(proj,ADCP.B10.Lat,ADCP.B10.Lon);
ADCP.B10.utm = [xutm,yutm];
ADCP.B10.depth = A.depth;

load(ADCPnames{2});
ADCP.B13.See = A.See;
ADCP.B13.Lat = A.latitude;
ADCP.B13.Lon = A.longitude;
[xutm,yutm] = projfwd(proj,ADCP.B13.Lat,ADCP.B13.Lon);
ADCP.B13.utm = [xutm,yutm];
ADCP.B13.depth = A.depth;

load(ADCPnames{3});
ADCP.X05.See = A.See;
ADCP.X05.Lat = A.latitude;
ADCP.X05.Lon = A.longitude;
[xutm,yutm] = projfwd(proj,ADCP.X05.Lat,ADCP.X05.Lon);
ADCP.X05.utm = [xutm,yutm];
ADCP.X05.depth = A.depth;

load(ADCPnames{4});
ADCP.X11.See = A.See;
ADCP.X11.Lat = A.latitude;
ADCP.X11.Lon = A.longitude;
[xutm,yutm] = projfwd(proj,ADCP.X11.Lat,ADCP.X11.Lon);
ADCP.X11.utm = [xutm,yutm];
ADCP.X11.depth = A.depth;


%% Plotting the time-averaged spectrum for the Smart Spotters

SMsee = {SM.E01.See, SM.E02.See, SM.E05.See, SM.E07.See, SM.E08.See, ...
    SM.E09.See, SM.E10.See, SM.E11.See, SM.E13.See};

figure(3);clf;
for i = 1:length(SMsee)
    subplot(3,1,1)
    SMAVGspec{i} = nanmean(SMsee{i},2);
    if i > 3
        subplot(3,1,2)
    end
    if i > 6
        subplot(3,1,3)
    end
    plot(SMf(1:200),SMAVGspec{i}(1:200),'LineWidth',1.5)
    hold on
    
end
subplot(3,1,1)
title('Time-Averaged Spectrum of the Smart Spotters')
legend('E01','E02','E05','location','best')
xlabel('Frequency (Hz)')
ylabel('Energy(m^2/Hz)')
ylim([0,0.7])
xlim([0,0.3])
grid on
subplot(3,1,2)
legend('E07','E08','E09','location','best')
grid on
xlabel('Frequency (Hz)')
ylabel('Energy (m^2/Hz)')
ylim([0,0.7])
xlim([0,0.3])
subplot(3,1,3)
legend('E10','E11','E13','location','best')
grid on
xlabel('Frequency (Hz)')
ylabel('Energy (m^2/Hz)')
ylim([0,0.7])
xlim([0,0.3])


SM.E01.AVGspec = SMAVGspec{1};
SM.E02.AVGspec = SMAVGspec{2};
SM.E05.AVGspec = SMAVGspec{3};
SM.E07.AVGspec = SMAVGspec{4};
SM.E08.AVGspec = SMAVGspec{5};
SM.E09.AVGspec = SMAVGspec{6};
SM.E10.AVGspec = SMAVGspec{7};
SM.E11.AVGspec = SMAVGspec{8};
SM.E13.AVGspec = SMAVGspec{9};


%% Plotting the time-averaged spectrum for the ADCPs (+ spotters)

ADCPsee = {ADCP.B10.See, ADCP.B13.See, ADCP.X05.See, ADCP.X11.See};


figure(4);clf;
subplot(2,1,1)
plot(Bfreq{1}(1:40),Bmeanspec{1}(1:40),'r')
hold on
plot(Bfreq{2}(1:40),Bmeanspec{2}(1:40),'b')
plot(Bfreq{3}(1:40),Bmeanspec{3}(1:40),'k')

subplot(2,1,2)
plot(Xfreq{1}(1:40),Xmeanspec{1}(1:40),'r')
hold on
plot(Xfreq{2}(1:40),Xmeanspec{2}(1:40),'b')
plot(Xfreq{3}(1:40),Xmeanspec{3}(1:40),'k')

for i = 1:length(ADCPsee)
    
    ADCPAVGspec{i} = nanmean(ADCPsee{i},2); %calculate the time averaged spectrum
    
    subplot(2,1,1) %China Rock plots go in this subplot
     if i > 2
         subplot(2,1,2) %Asilomar plots go in this subplot
     end
   plot(ADCPf(1:95),ADCPAVGspec{i}(1:95),'-','LineWidth',2) 
   hold on
end


subplot(2,1,1)
grid on
title('China Rock Time-Averaged Spectrum','fontsize',16)
xlabel('Frequency (Hz)','fontsize',13)
ylabel('Energy (m^2/Hz)','fontsize',13)
ylim([0 1.2])
legend('B01 (h=33.9m)','B03 (h=21.0m)','B05 (h=16.2m)',...
    'B10(ADCP) (h=9.5m)','B13(ADCP) (h=4.7m)','fontsize',10)

subplot(2,1,2)
grid on
title({'','Asilomar Time-Averaged Spectrum'},'fontsize',16)
xlabel('Frequency (Hz)','fontsize',13)
ylabel('Energy (m^2/Hz)','fontsize',13)
ylim([0 1.2])
legend('X01 (h=20.7m)','X03 (h=15.8m)','X04 (h=11.2m)',...
    'X05(ADCP) (h=9.3m)','X11(ADCP) (h=2.72m)','fontsize',10)


ADCP.B10.AVGspec = ADCPAVGspec{1};
ADCP.B13.AVGspec = ADCPAVGspec{2};
ADCP.X05.AVGspec = ADCPAVGspec{3};
ADCP.X11.AVGspec = ADCPAVGspec{4};



save('SM&ADCP_All.mat','SM','ADCP')



%% Creating Plots with Direction and Average Spectrums


    % China Rock
figure(5);clf;
subplot(2,1,1)
plot(Bfreq{1},BEMEMmean{1},'LineWidth',2)
hold on; grid on;
plot(Bfreq{2},BEMEMmean{2},'LineWidth',2)
plot(Bfreq{3},BEMEMmean{3},'LineWidth',2)
yline(BNormWaveDir) %add line for the shoreline normal direction
title('China Rock: Time-Averaged Directions')
xlabel('Freq (Hz)')
ylabel('Direction (degrees)')
legend('B01 (h = 33.9m)','B03 (h = 21.0m)','B05 (h = 16.2m)','Normal','location','southeast')
xlim([0 0.3])
ylim([225 325])

subplot(2,1,2)
plot(Bfreq{1}(1:40),Bmeanspec{1}(1:40),'LineWidth',2)
hold on; grid on;
plot(Bfreq{2}(1:40),Bmeanspec{2}(1:40),'LineWidth',2)
plot(Bfreq{3}(1:40),Bmeanspec{3}(1:40),'LineWidth',2)
plot(ADCPf(1:95),ADCP.B10.AVGspec(1:95),'LineWidth',2)
plot(ADCPf(1:95),ADCP.B13.AVGspec(1:95),'LineWidth',2)
title('China Rock: Time-Averaged Spectrums')
xlabel('Frequency (Hz)')
ylabel('Energy (m^2/Hz)')
xlim([0 0.3])
ylim([0 1.2])
legend('B01 (h = 33.9m)','B03 (h = 21.0m)','B05 (h = 16.2m)',...
    'B10 (h = 9.5m)','B13 (h = 4.7m)')


    % Asilomar
figure(6);clf;
subplot(2,1,1)
plot(Xfreq{1},XEMEMmean{1},'LineWidth',2)
hold on; grid on;
plot(Xfreq{2},XEMEMmean{2},'r','LineWidth',2)
plot(Xfreq{3},XEMEMmean{3},'LineWidth',2)
yline(XNormWaveDir) %add line for the shoreline normal direction
title('Asilomar: Time-Averaged Directions')
xlabel('Freq (Hz)')
ylabel('Direction (degrees)')
legend('X01 (h = 20.7m)','X03 (h = 15.8m)','X04 (h = 11.2m)','Normal','location','southeast')
xlim([0 0.3])
ylim([225 325])

subplot(2,1,2)
plot(Xfreq{1}(1:40),Xmeanspec{1}(1:40),'LineWidth',2)
hold on; grid on;
plot(Xfreq{2}(1:40),Xmeanspec{2}(1:40),'LineWidth',2)
plot(Xfreq{3}(1:40),Xmeanspec{3}(1:40),'LineWidth',2)
plot(ADCPf(1:95),ADCP.X05.AVGspec(1:95),'LineWidth',2)
plot(ADCPf(1:95),ADCP.X11.AVGspec(1:95),'LineWidth',2)
title('Asilomar: Time-Averaged Spectrums')
xlabel('Frequency (Hz)')
ylabel('Energy (m^2/Hz)')
xlim([0 0.3])
legend('X01 (h = 20.7m)','X03 (h = 15.8m)','X04 (h = 11.2m)',...
    'X05 (h = 9.3m)','X11 (h = 2.7m)')


%% B01 Avg Spec and avg (DEC PRES)



F = figure(6);clf;

subplot(2,1,1)
plot(Bfreq{3}(1:40),Bmeanspec{3}(1:40),'b','LineWidth',2.5);
hold on; grid on;
ax = gca(F);
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
xlabel('Frequency (Hz)','fontsize',22)
ylabel('Energy (m^2/Hz)','fontsize',22)
title('B05: Time-Averaged Energy Spectrum','fontsize',25)
ylim([0 1])
xlim([0 0.35])


subplot(2,1,2)
plot(Bfreq{3}(1:40),BEMEMmean{3}(1:40),'r','LineWidth',2.5)
hold on; grid on;
yline(BNormWaveDir,'linewidth',1.3)
title({[],['B05: Time-Averaged Wave Direction Spectrum']},'fontsize',25)
ax = gca(F);
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
xlabel('Frequency (Hz)','fontsize',22)
ylabel('Direction (^o)','fontsize',22)
ylim([240 310])
xlim([0 0.35])
legend('Spectrum','Shoreline Normal','fontsize',18,'location','SouthEast')



%% CR Buoy TA spectrums

figure(7);clf;
subplot(2,1,1)
plot(Bfreq{1}(1:40),Bmeanspec{1}(1:40),'r','LineWidth',2)
hold on
plot(Bfreq{2}(1:40),Bmeanspec{2}(1:40),'b','LineWidth',2)
plot(Bfreq{3}(1:40),Bmeanspec{3}(1:40),'g','LineWidth',2)

subplot(2,1,2)
plot(Xfreq{1}(1:40),Xmeanspec{1}(1:40),'r','LineWidth',2)
hold on
plot(Xfreq{2}(1:40),Xmeanspec{2}(1:40),'b','LineWidth',2)
plot(Xfreq{3}(1:40),Xmeanspec{3}(1:40),'g','LineWidth',2)


subplot(2,1,1)
grid on
title('China Rock Time-Averaged Spectrum','fontsize',16)
xlabel('Frequency (Hz)','fontsize',13)
ylabel('Energy (m^2/Hz)','fontsize',13)
xlim([0 0.35])
ylim([0 1.2])
legend('B01 (h=33.9m)','B03 (h=21.0m)','B05 (h=16.2m)','fontsize',10)

subplot(2,1,2)
grid on
title({'','Asilomar Time-Averaged Spectrum'},'fontsize',16)
xlabel('Frequency (Hz)','fontsize',13)
ylabel('Energy (m^2/Hz)','fontsize',13)
xlim([0 0.35])
ylim([0 1.2])
legend('X01 (h=20.7m)','X03 (h=15.8m)','X04 (h=11.2m)','fontsize',10)



