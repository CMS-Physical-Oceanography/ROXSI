%% plot_WBX01TS.m


    % Noah Clark
    % Created: 6/16/2023
    
    % Purpose: 
    %          - create a figure showing the time series plots of Hsig,
    %          energy-weighted mean and peak frequency, mean and peak
    %          bottom orbital velocity, mean and peak wavelength, and mean
    %          and peak wave orbital excursion for buoy X01 (plot #1)
    %          - calculate the mean and peak orbital excursion
    %          - create a 3 panel time-series figure of the spectrogram,
    %            wind speed, and tide data from NOAA for X01 (plot #2)
    %          - create the same plot as plot #2 but instead its datetime
    %            limits are from June 27th to July 6th          
    
    
%%

clear;clc;

    %Load in data from WBvariables.mat:
load('WBvariables.mat')


%%

                % Create plot #1 
x1 = datetime('15-Jun-2022 00:00:00','TimeZone','America/Los_Angeles');
x2 = datetime('21-Jul-2022 00:00:00','TimeZone','America/Los_Angeles');

figure(1);clf;
set(gcf,'position',[0,300,450,700])

    %Hsig:
subplot(5,1,1)
plot(Xtime{1},XGivenHsig{1},'k')
set(gca,'xticklabels',[])
title('X01 Time Series: Significant Wave Height')
ylabel({'H_s_i_g (m)',''},'fontsize',7.8)
ylim([0 2.5])
xlim([x1,x2])
grid on

    %Energy-Weighted Mean and Peak Frequency:
subplot(5,1,2)
plot(Xtime{1},fm)
hold on
plot(Xtime{1},fp,'r')
set(gca,'xticklabels',[])
title('X01 Time Series: Frequency')
ylabel({'Frequency (Hz)',''},'fontsize',7.8)
ylim([0 0.45])
xlim([x1,x2])
legend('Energy-Weighted Mean (f_m)','Peak (f_p)','Orientation',...
    'Horizontal','Location','north','fontsize',7)
grid on

    %Mean and Peak Bottom Orbital Velocity:
subplot(5,1,3)
plot(Xtime{1},mBOV)
hold on
plot(Xtime{1},pBOV,'r')
set(gca,'xticklabels',[])
title('X01 Time Series: Bottom Orbital Velocity')
ylabel({'Bottom Orbital Velocity (m/s)',''},'fontsize',7.8)
ylim([0 0.7])
xlim([x1,x2])
legend('Mean (u_m)','Peak (u_p)','Orientation','Horizontal',...
    'Location','north','fontsize',7)
grid on

    %Mean and Peak Wavelength:
subplot(5,1,4)
plot(Xtime{1},mWavelength)
hold on
plot(Xtime{1},pWavelength,'r')
set(gca,'xticklabels',[])
title('X01 Time Series: Wavelength')
ylabel({'Wavelength (m)',''},'fontsize',7.8)
ylim([0 375])
xlim([x1,x2])
legend('Mean (L_m)','Peak (L_p)','Orientation','Horizontal',...
    'Location','north','fontsize',7)
grid on


    % Calculate Wave Orbital Excursion:
%First must calculate the peak and mean (DOUBLE CHECK THAT I DID THIS RIGHT)
EXp = pBOV./(2*pi./Tp);
EXm = mBOV./(2*pi./Tm);
%Now, plotting
subplot(5,1,5)
plot(Xtime{1},EXm)
hold on
plot(Xtime{1},EXp,'r')
title('X01 Time Series: Wave Bottom Orbital Excursion')
ylabel({'Wave Orbital Excursion (m)',''},'fontsize',7.8)
xlabel('Date/Time')
ylim([0 1.65])
xlim([x1,x2])
legend('Mean (\xi_m)','Peak (\xi_p)','Orientation','Horizontal',...
    'Location','north','fontsize',7)
grid on



                % Create plot #2:
    %For the entire time of recording:
figure(2);clf;
set(gcf,'position',[450,300,500,700])

subplot(3,1,1)
pcolor(Xtime{1},Gfreq,XSee{1})
shading flat
cb = colorbar;
caxis([0 3])
set(cb,'position',[0.82 0.8429 0.0460 0.0658])
set(cb,'color','w')
ylim([0 0.4])
ylabel(cb,'Wave Energy (m^2/Hz)')
cb.Label.Color = 'k';
title('X01 Spectrogram');ylabel({'Frequency (Hz)',''},'fontsize',8);
xlim([x1,x2])
hold on

subplot(3,1,2)
plot(WindDT,WindSpd,'b')
hold on
plot(OffWindTime,OffWindSpd,'m')
title('Wind Speed in the Area');ylabel({'Wind Speed (m/s)',''},'fontsize',8);
xlim([x1,x2]);
grid on
xticks([datetime('15-Jun-2022 00:00:00','TimeZone','America/Los_Angeles'),...
    datetime('22-Jun-2022 00:00:00','TimeZone','America/Los_Angeles'),...
    datetime('29-Jun-2022 00:00:00','TimeZone','America/Los_Angeles'),...
    datetime('06-Jul-2022 00:00:00','TimeZone','America/Los_Angeles'),...
    datetime('13-Jul-2022 00:00:00','TimeZone','America/Los_Angeles'),...
    datetime('20-Jul-2022 00:00:00','TimeZone','America/Los_Angeles')])
xticklabels({'Jun 15','Jun 22','Jun 29','Jul 06','Jul 13','Jul 20'})
ylim([0 18]);
legend('Nearshore Wind','Offshore Wind','orientation','horizontal',...
    'location','north')

subplot(3,1,3)
plot(Xtime{1},NOAA.Verified_m_,'g')
title('Tide Data from NOAA Buoy (Station: 9413450)');
xlabel('Date');ylabel({'Water Elevation (m)(NAVD)',''},'fontsize',8);
xlim([x1,x2])
ylim([-1 2.5])
grid on


                    % Create plot #3:
    %From June 27th to July 6th:
x1_2 = datetime('27-Jun-2022 00:00:00','TimeZone','America/Los_Angeles');
x2_2 = datetime('06-Jul-2022 00:00:00','TimeZone','America/Los_Angeles');

figure(3);clf;
set(gcf,'position',[950,300,500,700])

subplot(3,1,1)
pcolor(Xtime{1},Gfreq,XSee{1})
shading flat
cb = colorbar;
caxis([0 3])
ylim([0 0.4])
set(cb,'position',[0.832 0.8429 0.0460 0.0658])
ylabel(cb,'Wave Energy (m^2/Hz)')
set(cb,'color','w')
cb.Label.Color = 'k';
title('X01 Spectrogram (Jun 27 - Jul 6)');ylabel({'Frequency (Hz)',''},'fontsize',8);
xlim([x1_2,x2_2])
hold on

subplot(3,1,2)
plot(WindDT,WindSpd,'b')
hold on
plot(OffWindTime,OffWindSpd,'m')
xlim([x1_2,x2_2])
ylim([0 16])
title('Wind Speed in the Area (Jun 27 - Jul 6)');ylabel({'Wind Speed (m/s)',''},'fontsize',8);
grid on
legend('Nearshore Wind','Offshore Wind','orientation','horizontal',...
    'location','north')

subplot(3,1,3)
plot(Xtime{1},NOAA.Verified_m_,'g')
title('Tide Data from NOAA Buoy (Station: 9413450) (Jun 27 - Jul 6)');
xlabel('Date');ylabel({'Water Elevation (m)(NAVD)',''},'fontsize',8);
xlim([x1_2,x2_2])
grid on



%%

    % Saving Variables:
save('WBvariables.mat','-append','EXp','EXm')





