%% plot_DailyShoalB.m
%
% Noah Clark        3/14/24
%
% Purpose: Calculate the observed and theoretical shoaling coefficients at
%          each hour and determine the EWM for each hour in the sea and 
%          swell bands. Then average the hours per day. This is done for
%          all spotter buoys and ADCPs at China Rock.
%
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


%% Preliminaries

clc;clear;

    % Spotter Buoy Data
load('WBvariables.mat','BSee','Bdepth','Btime')
h = Bdepth;
See = BSee;
See{1} = See{1}(2:end,:);
See{2} = See{2}(2:end,:);
See{3} = See{3}(2:end,:);

    % ADCP Data
load('SM&ADCP_All.mat','ADCP')
TdayADCP = datetime(2022,6,25):datetime(2022,7,19);
ADCP.time = ADCP.time(66:666); % June 25 - July 19
ADCP.B10.See = ADCP.B10.See(2:end,66:666); % take out the times to fit with data from spotters
ADCP.B13.See = ADCP.B13.See(2:end,66:666);
ADCP.B10.depth = ADCP.B10.depth(66:666);
ADCP.B13.depth = ADCP.B13.depth(66:666);

ADCP.freq = ADCP.freq(2:end);

    % Define Variables
ff = 0.0098:0.0098:(0.0098*128);
ff = ff';
TT = 1./ff;
ffa = ff(1:35); %for adcp
TTa = 1./ffa; %for adcp

Tday = datetime(2022,6,16):datetime(2022,7,19);
TdayB01 = [datetime(2022,6,16):datetime(2022,6,30),...
    datetime(2022,7,7):datetime(2022,7,19)];

    % Interpolate ADCP Sees to spotter frequencies
ADCP.B10.See = interp1(ADCP.freq',ADCP.B10.See,ff');
ADCP.B10.See = ADCP.B10.See(1:35,:);
ADCP.B13.See = interp1(ADCP.freq',ADCP.B13.See,ff');
ADCP.B13.See = ADCP.B13.See(1:35,:);


%% Calculate GROUP SPEED

% WAVELENGTH and PHASE SPEED
        % B01
for i = 1:725
    for j = 1:128
        [L{1}(j,i),~] = LDR(TT(j),h{1}(i)); %wavelength
    end

    c{1}(:,i) = L{1}(:,i)./TT; %phase speed

    for j = 1:128
        Cg{1}(j,i) = (c{1}(j,i)/2).*(1+((4*pi*h{1}(i))/L{1}(j,i))./...
            sinh((4*pi*h{1}(i))/L{1}(j,i)));
    end
end

        % B03 and B05
for i = 1:832
    for j = 1:128
        [L{2}(j,i),~] = LDR(TT(j),h{2}(i)); %wavelength
        [L{3}(j,i),~] = LDR(TT(j),h{3}(i)); %wavelength
    end

    c{2}(:,i) = L{2}(:,i)./TT; %phase speed
    c{3}(:,i) = L{3}(:,i)./TT; %phase speed

    for j = 1:128
        Cg{2}(j,i) = (c{2}(j,i)/2).*(1+((4*pi*h{2}(i))/L{2}(j,i))./...
            sinh((4*pi*h{2}(i))/L{2}(j,i)));
        Cg{3}(j,i) = (c{3}(j,i)/2).*(1+((4*pi*h{3}(i))/L{3}(j,i))./...
            sinh((4*pi*h{3}(i))/L{3}(j,i)));
    end
end

        % B10 and B13
for i = 1:601
    for j = 1:35
        [L{4}(j,i),~] = LDR(TTa(j),ADCP.B10.depth(i)); %wavelength
        [L{5}(j,i),~] = LDR(TTa(j),ADCP.B13.depth(i)); %wavelength
    end
    
    c{4}(:,i) = L{4}(:,i)./TTa; %phase speed
    c{5}(:,i) = L{5}(:,i)./TTa; %phase speed
    
    for j = 1:35
        Cg{4}(j,i) = (c{4}(j,i)/2).*(1+((4*pi*ADCP.B10.depth(i))/L{4}(j,i))./...
            sinh((4*pi*ADCP.B10.depth(i))/L{4}(j,i)));
        Cg{5}(j,i) = (c{5}(j,i)/2).*(1+((4*pi*ADCP.B13.depth(i))/L{5}(j,i))./...
            sinh((4*pi*ADCP.B13.depth(i))/L{5}(j,i)));
    end
end


%% Calculate Observed and Theoretical Shoaling Coeffs per hour

% THEORETICAL
Ks_T{1} = sqrt(Cg{1}./Cg{2}(:,[1:397 505:end]));
Ks_T{2} = sqrt(Cg{2}./Cg{3});
Ks_T{3} = sqrt(Cg{3}(1:35,217:817)./Cg{4});
Ks_T{4} = sqrt(Cg{4}./Cg{5});

% OBSERVED
Ks_O{1} = See{2}(:,[1:397 505:end])./See{1};
Ks_O{2} = See{3}./See{2};
Ks_O{3} = ADCP.B10.See./See{3}(1:35,217:817);
Ks_O{4} = ADCP.B13.See./ADCP.B10.See;

%% Determine EWM for each hour in Sea and Swell
% from here find the EWM for each hour in the sea and swell and then find
% the average of those per day 

% Average Sees between buoys
SeeB03Cut = See{2}(1:35,[1:397 505:end]); % trim out middle to match dates of B01 See
ASee{1} = (See{1}(1:35,:) + SeeB03Cut)./2;
ASee{2} = (See{2}(1:35,:) + See{3}(1:35,:))./2;
ASee{3} = (See{3}(1:35,217:817) + ADCP.B10.See(1:35,:))./2;
ASee{4} = (ADCP.B10.See(1:35,:) + ADCP.B13.See(1:35,:))./2;

for j = 1:4
    
    if j == 1 %B01 to B03
        N = 725;
    elseif j == 2 %B03 to B05
        N = 832;
    else %B05 to B10 & B10 to B13
        N = 601;
    end
    
    for i = 1:N
            % Sea
        m0Sea = trapz(ff(9:24),ASee{j}(9:24,i),1);
        m1SeaT = trapz(ff(9:24),ASee{j}(9:24,i).*Ks_T{j}(9:24,i));
        m1SeaO = trapz(ff(9:24),ASee{j}(9:24,i).*Ks_O{j}(9:24,i));
        EWM_KsT_Sea{j}(i) = m1SeaT./m0Sea;
        EWM_KsO_Sea{j}(i) = m1SeaO./m0Sea;
            % Swell
        m0Swell = trapz(ff(4:9),ASee{j}(4:9,i),1);
        m1SwellT = trapz(ff(4:9),ASee{j}(4:9,i).*Ks_T{j}(4:9,i));
        m1SwellO = trapz(ff(4:9),ASee{j}(4:9,i).*Ks_O{j}(4:9,i));
        EWM_KsT_Swell{j}(i) = m1SwellT./m0Swell;
        EWM_KsO_Swell{j}(i) = m1SwellO./m0Swell;
    end
end
        

%% Average per day

    % B01 to B03
% (**1st Section**)
for i = 1:15
    EWM_DD_KsT_Sea1(i) = mean(EWM_KsT_Sea{1}(i*24-23:i*24));
    EWM_DD_KsO_Sea1(i) = mean(EWM_KsO_Sea{1}(i*24-23:i*24));
    
    EWM_DD_KsT_Swell1(i) = mean(EWM_KsT_Swell{1}(i*24-23:i*24));
    EWM_DD_KsO_Swell1(i) = mean(EWM_KsO_Swell{1}(i*24-23:i*24));
end
% (**2nd Section**)
for i = 17.625:29.625
    EWM_DD_KsT_Sea2(i-16.625) = mean(EWM_KsT_Sea{1}(i*24-23:i*24));
    EWM_DD_KsO_Sea2(i-16.625) = mean(EWM_KsO_Sea{1}(i*24-23:i*24));
    
    EWM_DD_KsT_Swell2(i-16.625) = mean(EWM_KsT_Swell{1}(i*24-23:i*24));
    EWM_DD_KsO_Swell2(i-16.625) = mean(EWM_KsO_Swell{1}(i*24-23:i*24));
end

EWM_DD_KsO_Sea{1} = cat(2,EWM_DD_KsO_Sea1,EWM_DD_KsO_Sea2);
EWM_DD_KsO_Swell{1} = cat(2,EWM_DD_KsO_Swell1,EWM_DD_KsO_Swell2);
EWM_DD_KsT_Sea{1} = cat(2,EWM_DD_KsT_Sea1,EWM_DD_KsT_Sea2);
EWM_DD_KsT_Swell{1} = cat(2,EWM_DD_KsT_Swell1,EWM_DD_KsT_Swell2);

clear EWM_DD_KsO_Sea1 EWM_DD_KsO_Sea2 EWM_DD_KsO_Swell1 EWM_DD_KsO_Swell2
clear EWM_DD_KsT_Sea1 EWM_DD_KsT_Sea2 EWM_DD_KsT_Swell1 EWM_DD_KsT_Swell2


% B03 to B05, B05 to B10, & B10 to B13
for j = 2:4
    if j == 2 % B03 to B05
        N = 34;
    else % B05 to B10 & B10 to B13
        N = 25;
    end
    
    for i = 1:N
        EWM_DD_KsT_Sea{j}(i) = mean(EWM_KsT_Sea{j}(i*24-23:i*24));
        EWM_DD_KsO_Sea{j}(i) = mean(EWM_KsO_Sea{j}(i*24-23:i*24));
        
        EWM_DD_KsT_Swell{j}(i) = mean(EWM_KsT_Swell{j}(i*24-23:i*24));
        EWM_DD_KsO_Swell{j}(i) = mean(EWM_KsO_Swell{j}(i*24-23:i*24));
    end
end
    

%% Determine the means of the EWMs
% averages done for the observations

for i = 1:4
    % Sea
    Avg_KsO_Sea{i} = mean(EWM_DD_KsO_Sea{i});
    % Swell
    Avg_KsO_Swell{i} = mean(EWM_DD_KsO_Swell{i});
end


%% Plotting

figure(2);clf;

    % SWELL
subplot(1,2,1)
plot(TdayB01,EWM_DD_KsO_Swell{1},'bx','MarkerSize',10)
hold on; grid on;
plot(Tday,EWM_DD_KsO_Swell{2},'mx','MarkerSize',10)
plot(TdayADCP,EWM_DD_KsO_Swell{3},'gx','MarkerSize',10)
plot(TdayADCP,EWM_DD_KsO_Swell{4},'rx','MarkerSize',10)

plot(TdayB01,EWM_DD_KsT_Swell{1},'b-')
plot(Tday,EWM_DD_KsT_Swell{2},'m-')
plot(TdayADCP,EWM_DD_KsT_Swell{3},'g-')
plot(TdayADCP,EWM_DD_KsT_Swell{4},'r-')

title('China Rock: Swell Shoaling Coefficients','fontsize',17)
legend('B01 - B03 (AVG=1.534)','B03 - B05 (AVG=1.194)',...
    'B05 - B10 (AVG=0.953)','B10 - B13 (AVG=1.021)','fontsize',10,...
    'location','northwest')
xlabel('Date','fontsize',15)
ylabel('K_s','fontsize',15)
ylim([0.5 2.2])
xlim([datetime(2022,6,15) datetime(2022,7,20)])


    % SEA
subplot(1,2,2)
plot(TdayB01,EWM_DD_KsO_Sea{1},'bx','MarkerSize',10)
hold on; grid on;
plot(Tday,EWM_DD_KsO_Sea{2},'mx','MarkerSize',10)
plot(TdayADCP,EWM_DD_KsO_Sea{3},'gx','MarkerSize',10)
plot(TdayADCP,EWM_DD_KsO_Sea{4},'rx','MarkerSize',10)

plot(TdayB01,EWM_DD_KsT_Sea{1},'b-')
plot(Tday,EWM_DD_KsT_Sea{2},'m-')
plot(TdayADCP,EWM_DD_KsT_Sea{3},'g-')
plot(TdayADCP,EWM_DD_KsT_Sea{4},'r-')

title('China Rock: Sea Shoaling Coefficients','fontsize',17)
legend('B01 - B03 (AVG=0.943)','B03 - B05 (AVG=0.948)',...
    'B05 - B10 (AVG=1.113)','B10 - B13 (AVG=0.616)','fontsize',10,...
    'location','northwest')
xlabel('Date','fontsize',15)
ylim([0.5 2.2])
xlim([datetime(2022,6,15) datetime(2022,7,20)])



%% Plot Time Averaged for June 30th

figure(3);clf;
plot(ff,mean(See{3}(:,337:360),2),'-m')
hold on; grid on;
plot(ffa,mean(ADCP.B10.See(:,121:144),2),'-g')
plot(ffa,mean(ADCP.B13.See(:,121:144),2),'-r')

xline(0.0882,'--k')
xlim([0,0.25])
title('China Rock Time Averaged Spectrum for June 30th')
xlabel('f')
ylabel('E')
legend('B05','B10','B13')


%% Plot TIme Averaged for July  13th

figure(4);clf;
plot(ff,mean(See{3}(:,649:672),2),'-m')
hold on;grid on;
plot(ffa,mean(ADCP.B10.See(:,433:456),2),'-g')
plot(ffa,mean(ADCP.B13.See(:,433:456),2),'-r')

xline(0.0882,'--k')
xlim([0,0.25])
title('China Rock Time Averaged Spectrum for July 13th')
xlabel('f')
ylabel('E')
legend('B05','B10','B13')


%% Plot  Overall Time Averaged Spectrums

figure(5);clf;
plot(ff,mean(See{1},2),'b')
hold on; grid on;
plot(ff,nanmean(See{2},2),'m')
plot(ff,mean(See{3},2),'g')
plot(ffa,mean(ADCP.B10.See,2),'r')
plot(ffa,mean(ADCP.B13.See,2),'c')

load('SM&ADCP_All.mat','ADCP') %original ADCP Data
plot(ADCP.freq(2:end),nanmean(ADCP.B10.See(2:end,66:666),2),'r','LineWidth',2)
plot(ADCP.freq(2:end),nanmean(ADCP.B13.See(2:end,66:666),2),'c','LineWidth',2)

xline(0.0882,'--k')
xlim([0,0.25])
title('China Rock Time Averaged Spectrums')
xlabel('f')
ylabel('E')
legend('B01','B03','B05','B10','B13','og B10','og B13')



%% Clean Up

clear DD_Ks_O_1st DD_Ks_O_2nd DD_Ks_O_1st DD_Ks_t_2nd i j 
clear m0Sea m0Swell m1SeaO m1Seat m1SwellO m1SwellT N
















