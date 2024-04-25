%% plot_DailySumSeaSwellX.m
%
% Noah Clark            2/2/2024
%
%
% Purpose: Create two figures with 3 subplots each. One of the plots is the
%          total energy, the average direction, and the average directional
%          spread in the sea band per day for each of the 3 buoys at 
%          ASILOMAR. The other figure has the same 3 subplots but for the
%          swell band. 
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

%% General Preliminaries

clc;clear;

load('WBvariables.mat','XSee','Xtime','XEMEM','XNormWaveDir','Xdepth')

ff = 0:0.0098:(0.0098*128);
df = 0.0098;

addpath('../Summer (0516-0728)/Start Data')

% Time Vector (for each day)
Tday = datetime(2022,6,16):datetime(2022,7,19);


%% Calculate Snells Predicted


        % Sea Band - X01 to X03
for i = 1:832
    for j = 1:16
        % X01
        L1(j,i) = LDR(1/ff(j+9),Xdepth{1}(i));
        c1(j,i) = L1(j,i)/(1/ff(j+9));
        Cg1(j,i) = (c1(j,i)/2)*(1 + ((4*pi*Xdepth{1}(i)/L1(j,i))/sinh(4*pi*Xdepth{1}(i)/L1(j,i))));
        % X03
        L2(j,i) = LDR(1/ff(j+9),Xdepth{2}(i));
        c2(j,i) = L2(j,i)/(1/ff(j+9));
        Cg2(j,i) = (c2(j,i)/2)*(1 + ((4*pi*Xdepth{2}(i)/L2(j,i))/sinh(4*pi*Xdepth{2}(i)/L2(j,i))));
    end
end

Dir1 = XEMEM{1}(10:25,:) - XNormWaveDir;
DirPredSea{1} = asind(sind(Dir1).*(Cg2./Cg1)) + XNormWaveDir;
clear L1 c1 Cg1 L2 c2 Cg2 Dir1


        % Sea Band - X03 to X04
for i = 1:832
    for j = 1:16
        % X03
        L1(j,i) = LDR(1/ff(j+9),Xdepth{2}(i));
        c1(j,i) = L1(j,i)/(1/ff(j+9));
        Cg1(j,i) = (c1(j,i)/2)*(1 + ((4*pi*Xdepth{2}(i)/L1(j,i))/sinh(4*pi*Xdepth{2}(i)/L1(j,i))));
        % X04
        L2(j,i) = LDR(1/ff(j+9),Xdepth{3}(i));
        c2(j,i) = L2(j,i)/(1/ff(j+9));
        Cg2(j,i) = (c2(j,i)/2)*(1 + ((4*pi*Xdepth{3}(i)/L2(j,i))/sinh(4*pi*Xdepth{3}(i)/L2(j,i))));
    end
end

Dir1 = XEMEM{2}(10:25,:) - XNormWaveDir;
DirPredSea{2} = asind(sind(Dir1).*(Cg2./Cg1)) + XNormWaveDir;
clear L1 c1 Cg1 L2 c2 Cg2 Dir2



        % Swell Band - X01 to X03
for i = 1:832
    for j = 1:6
        % X01
        L1(j,i) = LDR(1/ff(j+4),Xdepth{1}(i));
        c1(j,i) = L1(j,i)/(1/ff(j+4));
        Cg1(j,i) = (c1(j,i)/2)*(1 + ((4*pi*Xdepth{1}(i)/L1(j,i))/sinh(4*pi*Xdepth{1}(i)/L1(j,i))));
        % X03
        L2(j,i) = LDR(1/ff(j+4),Xdepth{2}(i));
        c2(j,i) = L2(j,i)/(1/ff(j+4));
        Cg2(j,i) = (c2(j,i)/2)*(1 + ((4*pi*Xdepth{2}(i)/L2(j,i))/sinh(4*pi*Xdepth{2}(i)/L2(j,i))));
    end
end

Dir1 = XEMEM{1}(5:10,:) - XNormWaveDir;
DirPredSwell{1} = asind(sind(Dir1).*(Cg2./Cg1)) + XNormWaveDir;
clear L1 c1 Cg1 L2 c2 Cg2 Dir1


        % Swell Band - X03 to X04
for i = 1:832
    for j = 1:6
        % X03
        L1(j,i) = LDR(1/ff(j+4),Xdepth{2}(i));
        c1(j,i) = L1(j,i)/(1/ff(j+4));
        Cg1(j,i) = (c1(j,i)/2)*(1 + ((4*pi*Xdepth{2}(i)/L1(j,i))/sinh(4*pi*Xdepth{2}(i)/L1(j,i))));
        % X04
        L2(j,i) = LDR(1/ff(j+4),Xdepth{3}(i));
        c2(j,i) = L2(j,i)/(1/ff(j+4));
        Cg2(j,i) = (c2(j,i)/2)*(1 + ((4*pi*Xdepth{3}(i)/L2(j,i))/sinh(4*pi*Xdepth{3}(i)/L2(j,i))));
    end
end

Dir1 = XEMEM{2}(5:10,:) - XNormWaveDir;
DirPredSwell{2} = asind(sind(Dir1).*(Cg2./Cg1)) + XNormWaveDir;
clear L1 c1 Cg1 L2 c2 Cg2 Dir2


%% Calculate General


    % See    
SeaSee{1} = XSee{1}(10:25,:);
SwellSee{1} = XSee{1}(5:10,:);
SeaSee{2} = XSee{2}(10:25,:);
SwellSee{2} = XSee{2}(5:10,:);
SeaSee{3} = XSee{3}(10:25,:);
SwellSee{3} = XSee{3}(5:10,:);
% Sum all frequencies per hour
fSumSeaSee{1} = sum(SeaSee{1},1);
fSumSwellSee{1} = sum(SwellSee{1},1);
fSumSeaSee{2} = sum(SeaSee{2},1);
fSumSwellSee{2} = sum(SwellSee{2},1);
fSumSeaSee{3} = sum(SeaSee{3},1);
fSumSwellSee{3} = sum(SwellSee{3},1);

    % Directional Spread
load('roxsi_spotter_L2_X01_1151_reduced.mat')
Dspread{1} = spotterL2.EMEM.meandirspread_f;
load('roxsi_spotter_L2_X03_1157_reduced.mat')
Dspread{2} = spotterL2.EMEM.meandirspread_f;
load('roxsi_spotter_L2_X04_1155_reduced.mat')
Dspread{3} = spotterL2.EMEM.meandirspread_f;
clear spotterL2

    % Average/Sum Daily (everything)
for i = 1:34  % full 34 days:  16 June 2022 - 19 July 2022
    % Integrate See
    dSumSeaSee{1}(i) = nansum(fSumSeaSee{1}(i*24-23:i*24)).*df;
    dSumSwellSee{1}(i) = nansum(fSumSwellSee{1}(i*24-23:i*24)).*df;
    dSumSeaSee{2}(i) = nansum(fSumSeaSee{2}(i*24-23:i*24)).*df;
    dSumSwellSee{2}(i) = nansum(fSumSwellSee{2}(i*24-23:i*24)).*df;
    dSumSeaSee{3}(i) = nansum(fSumSeaSee{3}(i*24-23:i*24)).*df;
    dSumSwellSee{3}(i) = nansum(fSumSwellSee{3}(i*24-23:i*24)).*df;

        % EWM Direction and Predicted Direction for X01
    for j = i*24-23:i*24
        m0Sea = trapz(ff(10:25),XSee{1}(10:25,j),1);
        m1DirSea = trapz(ff(10:25),XSee{1}(10:25,j).*XEMEM{1}(10:25,j));
        m0Swell = trapz(ff(5:10),XSee{1}(5:10,j),1);
        m1DirSwell = trapz(ff(5:10),XSee{1}(5:10,j).*XEMEM{1}(5:10,j));
        m1DirSeaP = trapz(ff(10:25),XSee{1}(10:25,j).*DirPredSea{1}(:,j));
        m1DirSwellP = trapz(ff(5:10),XSee{1}(5:10,j).*DirPredSwell{1}(:,j));
        m1DirSpreadSea = trapz(ff(10:25),XSee{1}(10:25,j).*Dspread{1}(10:25,j));
        m1DirSpreadSwell = trapz(ff(5:10),XSee{1}(5:10,j).*Dspread{1}(5:10,j));
        
        EWM_dirSea(j-(i*24-24)) = m1DirSea./m0Sea;
        EWM_dirSwell(j-(i*24-24)) = m1DirSwell./m0Swell;
        EWM_dirSeaP(j-(i*24-24)) = m1DirSeaP./m0Sea;
        EWM_dirSwellP(j-(i*24-24)) = m1DirSwellP./m0Swell;
        EWM_dirSpreadSea(j-(i*24-24)) = m1DirSpreadSea./m0Sea;
        EWM_dirSpreadSwell(j-(i*24-24)) = m1DirSpreadSwell./m0Swell;
    end
    ADirSea{1}(i) = meanangle(EWM_dirSea) + 360;
    ADirSwell{1}(i) = meanangle(EWM_dirSwell) + 360;
    ADirPredSea{2}(i) = meanangle(EWM_dirSeaP) + 360;
    ADirPredSwell{2}(i) = meanangle(EWM_dirSwellP) + 360;
    ADspreadSea{1}(i) = meanangle(EWM_dirSpreadSea);
    ADspreadSwell{1}(i) = meanangle(EWM_dirSpreadSwell);


        % EWM Direction and Predicted Direction for X03
    for j = i*24-23:i*24
        m0Sea = trapz(ff(10:25),XSee{2}(10:25,j),1);
        m1DirSea = trapz(ff(10:25),XSee{2}(10:25,j).*XEMEM{2}(10:25,j));
        m0Swell = trapz(ff(5:10),XSee{2}(5:10,j),1);
        m1DirSwell = trapz(ff(5:10),XSee{2}(5:10,j).*XEMEM{2}(5:10,j));
        m1DirSeaP = trapz(ff(10:25),XSee{2}(10:25,j).*DirPredSea{2}(:,j));
        m1DirSwellP = trapz(ff(5:10),XSee{2}(5:10,j).*DirPredSwell{2}(:,j));
        m1DirSpreadSea = trapz(ff(10:25),XSee{2}(10:25,j).*Dspread{2}(10:25,j));
        m1DirSpreadSwell = trapz(ff(5:10),XSee{2}(5:10,j).*Dspread{2}(5:10,j));
        
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
    if i == 8 % to deal with the NaN at index 173 at X03
        ADirSea{2}(i) = meanangle(EWM_dirSea([1:4 6:end])) + 360;
        ADirSwell{2}(i) = meanangle(EWM_dirSwell([1:4 6:end])) + 360;
        ADirPredSea{3}(i) = meanangle(EWM_dirSeaP([1:4 6:end])) + 360;
        ADirPredSwell{3}(i) = meanangle(EWM_dirSwellP([1:4 6:end])) + 360;
        ADspreadSea{2}(i) = meanangle(EWM_dirSpreadSea([1:4 6:end]));
        ADspreadSwell{2}(i) = meanangle(EWM_dirSpreadSwell([1:4 6:end]));
    end
    

        % EWM Direction and Predicted Direction for X04
    for j = i*24-23:i*24
        m0Sea = trapz(ff(10:25),XSee{3}(10:25,j),1);
        m1DirSea = trapz(ff(10:25),XSee{3}(10:25,j).*XEMEM{3}(10:25,j));
        m0Swell = trapz(ff(5:10),XSee{3}(5:10,j),1);
        m1DirSwell = trapz(ff(5:10),XSee{3}(5:10,j).*XEMEM{3}(5:10,j));
        m1DirSpreadSea = trapz(ff(10:25),XSee{3}(10:25,j).*Dspread{3}(10:25,j));
        m1DirSpreadSwell = trapz(ff(5:10),XSee{3}(5:10,j).*Dspread{3}(5:10,j));
        
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
figure(3);clf;
subplot(3,1,1)
plot(Tday,dSumSeaSee{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,dSumSeaSee{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,dSumSeaSee{3},'xm','MarkerSize',10,'LineWidth',2)
ax = gca;
ax.FontSize = 15;
title('Sea Band: Integrated Energy','fontsize',23)
legend('X01','X03','X04','fontsize',14)
ylabel('Energy [m^2(???)]','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])

subplot(3,1,2)
plot(Tday,ADirSea{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,ADirSea{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,ADirSea{3},'xm','MarkerSize',10,'LineWidth',2)
plot(Tday,ADirPredSea{2},'-b','LineWidth',1)
plot(Tday,ADirPredSea{3},'-m','LineWidth',1)
yline(XNormWaveDir,'k')
ax = gca;
ax.FontSize = 15;
legend({'X01','X03','X04','X03 Predicted','X04 Predicted','Normal'},...
    'NumColumns',2,'fontsize',14,'location','southeast')
title('Sea Band: Average Direction','fontsize',23)
ylabel('Dir [^o]','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])
%ylim([235 310])

subplot(3,1,3)
plot(Tday,ADspreadSea{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,ADspreadSea{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,ADspreadSea{3},'xm','MarkerSize',10,'LineWidth',2)
ax = gca;
ax.FontSize = 15;
legend('X01','X03','X04','fontsize',14)
title('Sea Band: Average Directional Spread','fontsize',23)
ylabel('Spread [^o]','fontsize',21)
xlabel('Date','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])
%ylim([25 77])


    % Swell Plot:
figure(4);clf;
subplot(3,1,1)
plot(Tday,dSumSwellSee{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,dSumSwellSee{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,dSumSwellSee{3},'xm','MarkerSize',10,'LineWidth',2)
ax = gca;
ax.FontSize = 15;
legend('X01','X03','X04','fontsize',14)
title('Swell Band: Integrated Energy','fontsize',23)
ylabel('Energy [m^2(???)]','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])

subplot(3,1,2)
plot(Tday,ADirSwell{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,ADirSwell{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,ADirSwell{3},'xm','MarkerSize',10,'LineWidth',2)
plot(Tday,ADirPredSwell{2},'-b','LineWidth',1)
plot(Tday,ADirPredSwell{3},'-m','LineWidth',1)
yline(XNormWaveDir,'k')
ax = gca;
ax.FontSize = 15;
legend({'X01','X03','X04','X03 Predicted','X04 Predicted','Normal'},...
    'NumColumns',2,'fontsize',14)
title('Swell Band: Average Direction','fontsize',23)
ylabel('Dir [^o]','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])
%ylim([235 310])

subplot(3,1,3)
plot(Tday,ADspreadSwell{1},'xr','MarkerSize',10,'LineWidth',2)
hold on; grid on;
plot(Tday,ADspreadSwell{2},'xb','MarkerSize',10,'LineWidth',2)
plot(Tday,ADspreadSwell{3},'xm','MarkerSize',10,'LineWidth',2)
ax = gca;
ax.FontSize = 15;
legend('X01','X03','X04','fontsize',14)
title('Swell Band: Average Directional Spread','fontsize',23)
ylabel('Spread [^o]','fontsize',21)
xlabel('Date','fontsize',21)
xlim([datetime(2022,6,14) datetime(2022,7,21)])
%ylim([25 77])






