%% plot_DailySumSeaSwellB.m
%
% Noah Clark            2/2/2024
%
%
% Purpose: Create two figures with 3 subplots each. One of the plots is the
%          total energy, the average direction, and the average directional
%          spread in the sea band per day for each of the 3 buoys at CHINA
%          ROCK. The other figure has the same 3 subplots but for the swell
%          band. 
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

%% General Preliminaries

clc;clear;

load('WBvariables.mat','BSee','Btime','BEMEM','BNormWaveDir','Bdepth')

ff = 0:0.0098:(0.0098*128);
df = 0.0098;

addpath('../Summer (0516-0728)/Start Data')

% Time Vector (for each day) (used for B03 and B05)
Tday = datetime(2022,6,16):datetime(2022,7,19);
% Time Vector (for each day) (used only for B01)
%   - cut out middle section where buoy was not recording
TdayB01 = [datetime(2022,6,16):datetime(2022,6,30) datetime(2022,7,7):datetime(2022,7,19)];


%% Calculate Snells Predicted

        % Sea Band - B01 to B03
% B01
for i = 1:725
    for j = 1:16
        L1(j,i) = LDR(1/ff(j+9),Bdepth{1}(i));
        c1(j,i) = L1(j,i)/(1/ff(j+9));
        Cg1(j,i) = (c1(j,i)/2)*(1 + ((4*pi*Bdepth{1}(i)/L1(j,i))/sinh(4*pi*Bdepth{1}(i)/L1(j,i))));
    end
end

% B03 (to compare to B01)
B03D = Bdepth{2}([1:397 505:832]);
for i = 1:725
    for j = 1:16
        L2(j,i) = LDR(1/ff(j+9),B03D(i));
        c2(j,i) = L2(j,i)/(1/ff(j+9));
        Cg2(j,i) = (c2(j,i)/2)*(1 + ((4*pi*B03D(i)/L2(j,i))/sinh(4*pi*B03D(i)/L2(j,i))));
    end
end

Dir1 = BEMEM{1}(10:25,:) - BNormWaveDir;
DirPredSea{1} = asind(sind(Dir1).*(Cg2./Cg1)) + BNormWaveDir;
clear L1 c1 Cg1 L2 c2 Cg2 Dir1

        % Sea Band - B03 to B05
% B03
for i = 1:832
    for j = 1:16
        L1(j,i) = LDR(1/ff(j+9),Bdepth{2}(i));
        c1(j,i) = L1(j,i)/(1/ff(j+9));
        Cg1(j,i) = (c1(j,i)/2)*(1 + ((4*pi*Bdepth{2}(i)/L1(j,i))/sinh(4*pi*Bdepth{2}(i)/L1(j,i))));
    end
end

% B05
for i = 1:832
    for j = 1:16
        L2(j,i) = LDR(1/ff(j+9),Bdepth{3}(i));
        c2(j,i) = L2(j,i)/(1/ff(j+9));
        Cg2(j,i) = (c2(j,i)/2)*(1 + ((4*pi*Bdepth{3}(i)/L2(j,i))/sinh(4*pi*Bdepth{3}(i)/L2(j,i))));
    end
end

Dir1 = BEMEM{2}(10:25,:) - BNormWaveDir;
DirPredSea{2} = asind(sind(Dir1).*(Cg2./Cg1)) + BNormWaveDir;
clear L1 c1 Cg1 L2 c2 Cg2 Dir2




        % Swell Band - B01 to B03
% B01
for i = 1:725
    for j = 1:6
        L1(j,i) = LDR(1/ff(j+4),Bdepth{1}(i));
        c1(j,i) = L1(j,i)/(1/ff(j+4));
        Cg1(j,i) = (c1(j,i)/2)*(1 + ((4*pi*Bdepth{1}(i)/L1(j,i))/sinh(4*pi*Bdepth{1}(i)/L1(j,i))));
    end
end

% B03 (to compare to B01)
B03D = Bdepth{2}([1:397 505:832]);
for i = 1:725
    for j = 1:6
        L2(j,i) = LDR(1/ff(j+4),B03D(i));
        c2(j,i) = L2(j,i)/(1/ff(j+4));
        Cg2(j,i) = (c2(j,i)/2)*(1 + ((4*pi*B03D(i)/L2(j,i))/sinh(4*pi*B03D(i)/L2(j,i))));
    end
end

Dir1 = BEMEM{1}(5:10,:) - BNormWaveDir;
DirPredSwell{1} = asind(sind(Dir1).*(Cg2./Cg1)) + BNormWaveDir;
clear L1 c1 Cg1 L2 c2 Cg2 Dir1

        % Swell Band - B03 to B05
% B03
for i = 1:832
    for j = 1:6
        L1(j,i) = LDR(1/ff(j+4),Bdepth{2}(i));
        c1(j,i) = L1(j,i)/(1/ff(j+4));
        Cg1(j,i) = (c1(j,i)/2)*(1 + ((4*pi*Bdepth{2}(i)/L1(j,i))/sinh(4*pi*Bdepth{2}(i)/L1(j,i))));
    end
end

% B03 
for i = 1:832
    for j = 1:6
        L2(j,i) = LDR(1/ff(j+4),Bdepth{3}(i));
        c2(j,i) = L2(j,i)/(1/ff(j+4));
        Cg2(j,i) = (c2(j,i)/2)*(1 + ((4*pi*Bdepth{3}(i)/L2(j,i))/sinh(4*pi*Bdepth{3}(i)/L2(j,i))));
    end
end

Dir1 = BEMEM{2}(5:10,:) - BNormWaveDir;
DirPredSwell{2} = asind(sind(Dir1).*(Cg2./Cg1)) + BNormWaveDir;
clear L1 c1 Cg1 L2 c2 Cg2 Dir1


%% Calculate General
% - The frequency indecies for Sea: 0.088 - 0.245 Hz       (ff(10)-ff(25))
% - The frequency indecies for Swell: 0.039 - 0.088 Hz     (ff(5)-ff(10))


                            % B01

    % See
SeaSee{1} = BSee{1}(10:25,:);
SwellSee{1} = BSee{1}(5:10,:);
% Sum all frequencies per hour
fSumSeaSee{1} = sum(SeaSee{1},1);
fSumSwellSee{1} = sum(SwellSee{1},1);

    % Directional Spread
load('roxsi_spotter_L2_B01_1158_reduced.mat')
Dspread1 = spotterL2.EMEM.meandirspread_f;
load('roxsi_spotter_L2_B01_1150_reduced.mat')
Dspread2 = spotterL2.EMEM.meandirspread_f;
Dspread{1} = cat(2,Dspread1,Dspread2);

% Average/Sum Daily (everything) (** 1st Section **)
for i = 1:15 % first 15 days: 16 June 2022 - 30 June 2022
    % Integrate See
    dSumSeaSee1(i) = sum(fSumSeaSee{1}(i*24-23:i*24)).*df; %4*sqrt??? (H)
    dSumSwellSee1(i) = sum(fSumSwellSee{1}(i*24-23:i*24)).*df;

    % EWM Direction, Predicted Direction, and Directional Spread
    for j = i*24-23:i*24
        m0Sea = trapz(ff(10:25),BSee{1}(10:25,j),1);
        m1DirSea = trapz(ff(10:25),BSee{1}(10:25,j).*BEMEM{1}(10:25,j));
        m0Swell = trapz(ff(5:10),BSee{1}(5:10,j),1);
        m1DirSwell = trapz(ff(5:10),BSee{1}(5:10,j).*BEMEM{1}(5:10,j));
        m1DirSeaP = trapz(ff(10:25),BSee{1}(10:25,j).*DirPredSea{1}(:,j));
        m1DirSwellP = trapz(ff(5:10),BSee{1}(5:10,j).*DirPredSwell{1}(:,j));
        m1DirSpreadSea = trapz(ff(10:25),BSee{1}(10:25,j).*Dspread{1}(10:25,j));
        m1DirSpreadSwell = trapz(ff(5:10),BSee{1}(5:10,j).*Dspread{1}(5:10,j));
        
        EWM_dirSea(j-(i*24-24)) = m1DirSea./m0Sea;
        EWM_dirSwell(j-(i*24-24)) = m1DirSwell./m0Swell;
        EWM_dirSeaP(j-(i*24-24)) = m1DirSeaP./m0Sea;
        EWM_dirSwellP(j-(i*24-24)) = m1DirSwellP./m0Swell;
        EWM_dirSpreadSea(j-(i*24-24)) = m1DirSpreadSea./m0Sea;
        EWM_dirSpreadSwell(j-(i*24-24)) = m1DirSpreadSwell./m0Swell;
    end
    ADirSea1(i) = meanangle(EWM_dirSea) + 360;
    ADirSwell1(i) = meanangle(EWM_dirSwell) + 360;
    ADirPredSea1(i) = meanangle(EWM_dirSeaP) + 360;
    ADirPredSwell1(i) = meanangle(EWM_dirSwellP) + 360;
    ADirSpreadSea1(i) = meanangle(EWM_dirSpreadSea);
    ADirSpreadSwell1(i) = meanangle(EWM_dirSpreadSwell);

end

clear EWM_dirSea EWM_dirSwell EWM_dirSeaP EWM_dirSwellP EWM_dirSpreadSea EWM_dirSpreadSwell

% Average/Sum Daily (everything) (** 2nd Section **)
for i =  17.625:29.625 % last 13 days: 07 July 2022 - 19 July 2022
    % Integrate See
    dSumSeaSee2(i-16.625) = sum(fSumSeaSee{1}(i*24-23:i*24)).*df;
    dSumSwellSee2(i-16.625) = sum(fSumSwellSee{1}(i*24-23:i*24)).*df;

        % EWM Direction and Predicted Direction
    for j = i*24-23:i*24
        m0Sea = trapz(ff(10:25),BSee{1}(10:25,j),1);
        m1DirSea = trapz(ff(10:25),BSee{1}(10:25,j).*BEMEM{1}(10:25,j));
        m0Swell = trapz(ff(5:10),BSee{1}(5:10,j),1);
        m1DirSwell = trapz(ff(5:10),BSee{1}(5:10,j).*BEMEM{1}(5:10,j));
        m1DirSeaP = trapz(ff(10:25),BSee{1}(10:25,j).*DirPredSea{1}(:,j));
        m1DirSwellP = trapz(ff(5:10),BSee{1}(5:10,j).*DirPredSwell{1}(:,j));
        m1DirSpreadSea = trapz(ff(10:25),BSee{1}(10:25,j).*Dspread{1}(10:25,j));
        m1DirSpreadSwell = trapz(ff(5:10),BSee{1}(5:10,j).*Dspread{1}(5:10,j));
        
        EWM_dirSea(j-(i*24-24)) = m1DirSea./m0Sea;
        EWM_dirSwell(j-(i*24-24)) = m1DirSwell./m0Swell;
        EWM_dirSeaP(j-(i*24-24)) = m1DirSeaP./m0Sea;
        EWM_dirSwellP(j-(i*24-24)) = m1DirSwellP./m0Swell;
        EWM_dirSpreadSea(j-(i*24-24)) = m1DirSpreadSea./m0Sea;
        EWM_dirSpreadSwell(j-(i*24-24)) = m1DirSpreadSwell./m0Swell;
    end
    ADirSea2(i-16.625) = meanangle(EWM_dirSea) + 360;
    ADirSwell2(i-16.625) = meanangle(EWM_dirSwell) + 360;
    ADirPredSea2(i-16.625) = meanangle(EWM_dirSeaP) + 360;
    ADirPredSwell2(i-16.625) = meanangle(EWM_dirSwellP) + 360;
    ADirSpreadSea2(i-16.625) = meanangle(EWM_dirSpreadSea);
    ADirSpreadSwell2(i-16.625) = meanangle(EWM_dirSpreadSwell);

end

clear EWM_dirSea EWM_dirSwell EWM_dirSeaP EWM_dirSwellP EWM_dirSpreadSea EWM_dirSpreadSwell

% Concatonate sections
dSumSeaSee{1} = cat(2,dSumSeaSee1,dSumSeaSee2);
dSumSwellSee{1} = cat(2,dSumSwellSee1,dSumSwellSee2);
ADirSea{1} = cat(2,ADirSea1,ADirSea2);
ADirSwell{1} = cat(2,ADirSwell1,ADirSwell2);
ADspreadSea{1} = cat(2,ADirSpreadSea1,ADirSpreadSea2);
ADspreadSwell{1} = cat(2,ADirSpreadSwell1,ADirSpreadSwell2);
ADirPredSea{2} = cat(2,ADirPredSea1,ADirPredSea2);
ADirPredSwell{2} = cat(2,ADirPredSwell1,ADirPredSwell2);
clear dSumSeaSee1 dSumSeaSee2 dSumSwellSee1 dSumSwellSee2 ADirSea1 ...
    ADirSea2 ADirSwell1 ADirSwell2 ADirSpreadSea1 ADirSpreadSea2...
    ADirSpreadSwell1 ADirSpreadSwell2 m0Sea m0Swell


%%
                            % B03 and B05

    % See    
SeaSee{2} = BSee{2}(10:25,:);
SwellSee{2} = BSee{2}(5:10,:);
SeaSee{3} = BSee{3}(10:25,:);
SwellSee{3} = BSee{3}(5:10,:);
% Sum all frequencies per hour
fSumSeaSee{2} = sum(SeaSee{2},1);
fSumSwellSee{2} = sum(SwellSee{2},1);
fSumSeaSee{3} = sum(SeaSee{3},1);
fSumSwellSee{3} = sum(SwellSee{3},1);

    % Directional Spread
load('roxsi_spotter_L2_B03_1152_reduced.mat')
Dspread{2} = spotterL2.EMEM.meandirspread_f;
load('roxsi_spotter_L2_B05_1153_reduced.mat')
Dspread{3} = spotterL2.EMEM.meandirspread_f;
clear spotterL2

    % Average/Sum Daily (everything)
for i = 1:34  % full 34 days:  16 June 2022 - 19 July 2022
    % Integrate See
    dSumSeaSee{2}(i) = sum(fSumSeaSee{2}(i*24-23:i*24)).*df;
    dSumSwellSee{2}(i) = sum(fSumSwellSee{2}(i*24-23:i*24)).*df;
    dSumSeaSee{3}(i) = sum(fSumSeaSee{3}(i*24-23:i*24)).*df;
    dSumSwellSee{3}(i) = sum(fSumSwellSee{3}(i*24-23:i*24)).*df;



        % EWM Direction and Predicted Direction for B03
    for j = i*24-23:i*24
        m0Sea = trapz(ff(10:25),BSee{2}(10:25,j),1);
        m1DirSea = trapz(ff(10:25),BSee{2}(10:25,j).*BEMEM{2}(10:25,j));
        m0Swell = trapz(ff(5:10),BSee{2}(5:10,j),1);
        m1DirSwell = trapz(ff(5:10),BSee{2}(5:10,j).*BEMEM{2}(5:10,j));
        m1DirSeaP = trapz(ff(10:25),BSee{2}(10:25,j).*DirPredSea{2}(:,j));
        m1DirSwellP = trapz(ff(5:10),BSee{2}(5:10,j).*DirPredSwell{2}(:,j));
        m1DirSpreadSea = trapz(ff(10:25),BSee{2}(10:25,j).*Dspread{2}(10:25,j));
        m1DirSpreadSwell = trapz(ff(5:10),BSee{2}(5:10,j).*Dspread{2}(5:10,j));
        
        EWM_dirSea(j-(i*24-24)) = m1DirSea./m0Sea;
        EWM_dirSwell(j-(i*24-24)) = m1DirSwell./m0Swell;
        EWM_dirSeaP(j-(i*24-24)) = m1DirSeaP./m0Sea;
        EWM_dirSwellP(j-(i*24-24)) = m1DirSwellP./m0Swell;
        EWM_dirSpreadSea(j-(i*24-24)) = m1DirSpreadSea./m0Sea;
        EWM_dirSpreadSwell(j-(i*24-24)) = m1DirSpreadSwell./m0Swell;
    end
    ADirSea{2}(i) = meanangle(EWM_dirSea) + 360;
    ADirSwell{2}(i) = meanangle(EWM_dirSwell) + 360;
    ADirPredSea{3}(i) = meanangle(EWM_dirSeaP) + 360;
    ADirPredSwell{3}(i) = meanangle(EWM_dirSwellP) + 360;
    ADspreadSea{2}(i) = meanangle(EWM_dirSpreadSea);
    ADspreadSwell{2}(i) = meanangle(EWM_dirSpreadSwell);


            % EWM Direction and Predicted Direction for B05
    for j = i*24-23:i*24
        m0Sea = trapz(ff(10:25),BSee{3}(10:25,j),1);
        m1DirSea = trapz(ff(10:25),BSee{3}(10:25,j).*BEMEM{3}(10:25,j));
        m0Swell = trapz(ff(5:10),BSee{3}(5:10,j),1);
        m1DirSwell = trapz(ff(5:10),BSee{3}(5:10,j).*BEMEM{3}(5:10,j));
        m1DirSpreadSea = trapz(ff(10:25),BSee{3}(10:25,j).*Dspread{3}(10:25,j));
        m1DirSpreadSwell = trapz(ff(5:10),BSee{3}(5:10,j).*Dspread{3}(5:10,j));
        
        EWM_dirSea(j-(i*24-24)) = m1DirSea./m0Sea;
        EWM_dirSwell(j-(i*24-24)) = m1DirSwell./m0Swell;
        EWM_dirSpreadSea(j-(i*24-24)) = m1DirSpreadSea./m0Sea;
        EWM_dirSpreadSwell(j-(i*24-24)) = m1DirSpreadSwell./m0Swell;
    end
    ADirSea{3}(i) = meanangle(EWM_dirSea) + 360;
    ADirSwell{3}(i) = meanangle(EWM_dirSwell) + 360;
    ADspreadSea{3}(i) = meanangle(EWM_dirSpreadSea);
    ADspreadSwell{3}(i) = meanangle(EWM_dirSpreadSwell);

end


%% Plot

    % Sea Plot:
figure(1);clf;
subplot(3,1,1)
plot(TdayB01,dSumSeaSee{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,dSumSeaSee{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,dSumSeaSee{3},'xm','MarkerSize',10,'LineWidth',2)
ax = gca;
ax.FontSize = 15;
title('Sea Band: Integrated Energy','fontsize',23)
legend('B01','B03','B05','fontsize',14)
ylabel('Energy [m^2(???)]','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])

subplot(3,1,2)
plot(TdayB01,ADirSea{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,ADirSea{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,ADirSea{3},'xm','MarkerSize',10,'LineWidth',2)
plot(TdayB01(1:15),ADirPredSea{2}(1:15),'-b','LineWidth',1)
plot(TdayB01(16:end),ADirPredSea{2}(16:end),'-b','LineWidth',1)
plot(Tday,ADirPredSea{3},'-m','LineWidth',1)
yline(BNormWaveDir,'k')
ax = gca;
ax.FontSize = 15;
legend({'B01','B03','B05','B03 Predicted','','B05 Predicted','Normal'},...
    'NumColumns',2,'fontsize',14,'location','southeast')
title('Sea Band: Average Direction','fontsize',23)
ylabel('Dir [^o]','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])
ylim([235 310])

subplot(3,1,3)
plot(TdayB01,ADspreadSea{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,ADspreadSea{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,ADspreadSea{3},'xm','MarkerSize',10,'LineWidth',2)
ax = gca;
ax.FontSize = 15;
legend('B01','B03','B05','fontsize',14)
title('Sea Band: Average Directional Spread','fontsize',23)
ylabel('Spread [^o]','fontsize',21)
xlabel('Date','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])
ylim([25 77])


    % Swell Plot:
figure(2);clf;
subplot(3,1,1)
plot(TdayB01,dSumSwellSee{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,dSumSwellSee{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,dSumSwellSee{3},'xm','MarkerSize',10,'LineWidth',2)
ax = gca;
ax.FontSize = 15;
legend('B01','B03','B05','fontsize',14)
title('Swell Band: Integrated Energy','fontsize',23)
ylabel('Energy [m^2(???)]','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])

subplot(3,1,2)
plot(TdayB01,ADirSwell{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,ADirSwell{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,ADirSwell{3},'xm','MarkerSize',10,'LineWidth',2)
plot(TdayB01(1:15),ADirPredSwell{2}(1:15),'-b','LineWidth',1)
plot(TdayB01(16:end),ADirPredSwell{2}(16:end),'-b','LineWidth',1)
plot(Tday,ADirPredSwell{3},'-m','LineWidth',1)
yline(BNormWaveDir,'k')
ax = gca;
ax.FontSize = 15;
legend({'B01','B03','B05','B03 Predicted','','B05 Predicted','Normal'},...
    'NumColumns',2,'fontsize',14)
title('Swell Band: Average Direction','fontsize',23)
ylabel('Dir [^o]','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])
ylim([235 310])

subplot(3,1,3)
plot(TdayB01,ADspreadSwell{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,ADspreadSwell{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,ADspreadSwell{3},'xm','MarkerSize',10,'LineWidth',2)
ax = gca;
ax.FontSize = 15;
legend('B01','B03','B05','fontsize',14)
title('Swell Band: Average Directional Spread','fontsize',23)
ylabel('Spread [^o]','fontsize',21)
xlabel('Date','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])
ylim([25 77])


%% Save Figures

% Change the directory to where I want to save the figures
% cd 'Spring Figures/Daily_Averages'
% saveas(figure(1),'Sea_DailyAvgs.jpeg')
% saveas(figure(2),'Swell_DailyAvgs.jpeg')








