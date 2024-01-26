%% calc_Snells.m
%
% Noah Clark        11/10/23
%
%
% Purpose: Use Snell's Law to determine the angle of refraction at
%           different frequency bands in the peak frequency range
%
%


%% Preliminaries

clc;clear;

load('WBvariables.mat','XSee','BSee','XNormWaveDir','BNormWaveDir',...
    'XEMEM','BEMEM','Bmeanspec','Xmeanspec','Xfreq','Bfreq',...
    'Xdepth','Bdepth')
load('SM&ADCP_All.mat') 

ff = 0.0098:0.0098:0.9961; %Hz
g = 9.81; %m/s^2
df = 0.0098; %Hz

for i = 1:3
    H_SR{i} = mean(4*sqrt(BSee{i}(7:9,:).*df));
    H_LR{i} = mean(4*sqrt(BSee{i}(11:15,:).*df));
end


%% Calculate the different Frequency Bands
% 75% of peak energy

    % 1st (smaller) frequency band:
[m1,indm1] = max(Bmeanspec{3}(1:9));
BoundE1 = m1 - (0.25*m1); % = 0.3004
indminf1 = 7; % f = 0.0686 Hz       E = 0.3412
indmaxf1 = 9; % f = 0.0882 Hz       E = 0.3065


    % 2nd (bigger) frequency band:
[m2,indm2] = max(Bmeanspec{3});
BoundE2 = m2 - (0.25*m2); % = 0.6220
indminf2 = 11; % f = 0.1078 Hz      E = 0.6086
indmaxf2 = 15; % f = 0.1470 Hz      E = 0.6287


%% From B01 to B03 (Time Averaged)
%         
% 
% % Average Depths:        
% BmeanDepth1 = mean(Bdepth{1});
% BmeanDepth2 = mean(Bdepth{2});
% 
% % Average Direction (over frequency):
% Bdirec1 = meanangle(BEMEM{1},2) + 360 - 180; %direction they are going to
% Bdirec2 = meanangle(BEMEM{2},2) + 360 - 180;
% 
% BNormDir = BNormWaveDir - 180; %direction they are going to
% 
% ObsTheta1 = BNormDir - Bdirec1; %wave directions in reference to shoreline
% ObsTheta2 = BNormDir - Bdirec2;
% 
% % Calculate Energy Weighted Mean Direction:
% EWdir_SR = sum(ObsTheta1(7:9)'.*Bmeanspec{1}(7:9))/sum(Bmeanspec{1}(7:9)); %for small freq range
% EWdir_LR = sum(ObsTheta1(11:15)'.*Bmeanspec{1}(11:15))/sum(Bmeanspec{1}(11:15)); %for large freq range
% 
% 
%     % 1st Frequency Band:
% [L1,~] = LDR(1/ff(8),BmeanDepth1); %use ff(8) because it is the middle of the 3 freqs
% c = L1/(1/ff(8));
% Cg1 = (c/2)*(1 + ((4*pi*BmeanDepth1/L1)/sinh(4*pi*BmeanDepth1/L1)));
% 
% [L2,~] = LDR(1/ff(8),BmeanDepth2); %use ff(8) because it is the middle of the 3 freqs
% c = L2/(1/ff(8));
% Cg2 = (c/2)*(1 + ((4*pi*BmeanDepth2/L2)/sinh(4*pi*BmeanDepth2/L2)));
% 
% %Snell:
% theta2 = asind((Cg2/Cg1)*sind(EWdir_SR))
% Compare = meanangle(ObsTheta2(7:9))
% 
% % 2nd Frequency Band
% [L1,~] = LDR(1/ff(13),BmeanDepth1); %use ff(8) because it is the middle of the 3 freqs
% c = L1/(1/ff(13));
% Cg1 = (c/2)*(1 + ((4*pi*BmeanDepth1/L1)/sinh(4*pi*BmeanDepth1/L1)));
% 
% [L2,~] = LDR(1/ff(13),BmeanDepth2); %use ff(8) because it is the middle of the 3 freqs
% c = L2/(1/ff(13));
% Cg2 = (c/2)*(1 + ((4*pi*BmeanDepth2/L2)/sinh(4*pi*BmeanDepth2/L2)));
% 
% %Snell:
% theta2 = asind((Cg2/Cg1)*sind(EWdir_LR))
% Compare = meanangle(ObsTheta2(11:15))




%% From B01 to B03 (Time Varying)
        
% Assign Energy Spectra:
See1 = BSee{1};
See2 = BSee{2}(:,[1:397 505:end]);

% Average Depths:        
Depth1 = Bdepth{1};
Depth2 = Bdepth{2}([1:397 505:end]);

% Direction (over frequency):
Bdirec1_SR = BEMEM{1}(7:9,:) - 180; %direction they are going to
Bdirec1_LR = BEMEM{1}(11:15,:) - 180;
Bdirec2_SR = BEMEM{2}(7:9,[1:397 505:end]) - 180;
Bdirec2_LR = BEMEM{2}(11:15,[1:397 505:end]) - 180;

BNormDir = BNormWaveDir - 180; %direction they are going to

ObsTheta1_SR = BNormDir - Bdirec1_SR; %wave directions in reference to shoreline
ObsTheta1_LR = BNormDir - Bdirec1_LR;
ObsTheta2_SR = BNormDir - Bdirec2_SR;
ObsTheta2_LR = BNormDir - Bdirec2_LR;

% Calculate Energy Weighted Mean Direction:
for i = 1:length(See1) 
    EWdir_SR1(i) = sum(ObsTheta1_SR(:,i).*See1(7:9,i))./sum(See1(7:9,i)); %for small freq range
    EWdir_LR1(i) = sum(ObsTheta1_LR(:,i).*See1(11:15,i))./sum(See1(11:15,i)); %for large freq range

    EWdir_SR2(i) = sum(ObsTheta2_SR(:,i).*See2(7:9,i))./sum(See2(7:9,i)); %for small freq range
    EWdir_LR2(i) = sum(ObsTheta2_LR(:,i).*See2(11:15,i))./sum(See2(11:15,i)); %for large freq range
end

    % 1st Frequency Band:
for i = 1:length(See1)
    [L1(i),~] = LDR(1/ff(8),Depth1(i)); %use ff(8) because it is the middle of the 3 freqs
    c1(i) = L1(i)/(1/ff(8));
    Cg1(i) = (c1(i)/2)*(1 + ((4*pi*Depth1(i)/L1(i))/sinh(4*pi*Depth1(i)/L1(i))));

    [L2(i),~] = LDR(1/ff(8),Depth2(i)); %use ff(8) because it is the middle of the 3 freqs
    c2(i) = L2(i)/(1/ff(8));
    Cg2(i) = (c2(i)/2)*(1 + ((4*pi*Depth2(i)/L2(i))/sinh(4*pi*Depth2(i)/L2(i))));
end


    % Calculate the wave height using shoaling and ref coeff:
H_SR_snell{1} = H_SR{1}./(sqrt(Cg2./Cg1).*sqrt(cosd(EWdir_SR1)./cosd(EWdir_SR2))); %snell height at B03 (SR)


% Snell:
B01obs_SR = meanangle(ObsTheta1_SR);
theta2_1stBand{1} = asind((Cg2./Cg1).*sind(EWdir_SR1));
Compare_1stBand{1} = meanangle(ObsTheta2_SR);
% Plot:
figure(1);clf; 
plot(theta2_1stBand{1},'r','LineWidth',2); hold on; grid on;
plot(Compare_1stBand{1},'b','LineWidth',2);
title('B03: Snells Law vs Observed (Freq: 0.0686-0.0882 Hz)')
legend('Snells Dir','Obs Dir')
ylabel('Wave Dir (relative to shore normal)')



    % 2nd Frequency Band:
for i = 1:length(See1)
    [L1(i),~] = LDR(1/ff(13),Depth1(i)); %use ff(13) because it is the middle of the 5 freqs
    c1(i) = L1(i)/(1/ff(13));
    Cg1(i) = (c1(i)/2)*(1 + ((4*pi*Depth1(i)/L1(i))/sinh(4*pi*Depth1(i)/L1(i))));

    [L2(i),~] = LDR(1/ff(13),Depth2(i)); %use ff(13) because it is the middle of the 5 freqs
    c2(i) = L2(i)/(1/ff(13));
    Cg2(i) = (c2(i)/2)*(1 + ((4*pi*Depth2(i)/L2(i))/sinh(4*pi*Depth2(i)/L2(i))));
end


    % Calculate the wave height using shoaling and ref coeff:
H_LR_snell{1} = H_LR{1}./(sqrt(Cg2./Cg1).*sqrt(cosd(EWdir_LR1)./cosd(EWdir_LR2))); %snell height at B03 (LR)


% Snell:
theta2_2ndBand{1} = asind((Cg2./Cg1).*sind(EWdir_LR1));
Compare_2ndBand{1} = meanangle(ObsTheta2_LR);
% Plot:
figure(2);clf; 
plot(theta2_2ndBand{1},'r','LineWidth',2); hold on; grid on;
plot(Compare_2ndBand{1},'b','LineWidth',2);
title('B03: Snells Law vs Observed (Freq: 0.1078-0.1470 Hz)')
legend('Snells Dir','Obs Dir')
ylabel('Wave Dir (relative to shore normal)')



%% B03 to B05

   
% Assign Energy Spectra:
See1 = BSee{2};
See2 = BSee{3};

% Average Depths:        
Depth1 = Bdepth{2};
Depth2 = Bdepth{3};

% Direction (over frequency):
Bdirec1_SR = BEMEM{2}(7:9,:) - 180; %direction they are going to
Bdirec1_LR = BEMEM{2}(11:15,:) - 180;
Bdirec2_SR = BEMEM{3}(7:9,:) - 180;
Bdirec2_LR = BEMEM{3}(11:15,:) - 180;

BNormDir = BNormWaveDir - 180; %direction they are going to

ObsTheta1_SR = BNormDir - Bdirec1_SR; %wave directions in reference to shoreline
ObsTheta1_LR = BNormDir - Bdirec1_LR;
ObsTheta2_SR = BNormDir - Bdirec2_SR;
ObsTheta2_LR = BNormDir - Bdirec2_LR;

% Calculate Energy Weighted Mean Direction:
for i = 1:length(See1) 
    EWdir_SR1(i) = sum(ObsTheta1_SR(:,i).*See1(7:9,i))./sum(See1(7:9,i)); %for small freq range
    EWdir_LR1(i) = sum(ObsTheta1_LR(:,i).*See1(11:15,i))./sum(See1(11:15,i)); %for large freq range
    EWdir_SR2(i) = sum(ObsTheta2_SR(:,i).*See2(7:9,i))./sum(See2(7:9,i)); %for small freq range
    EWdir_LR2(i) = sum(ObsTheta2_LR(:,i).*See2(11:15,i))./sum(See2(11:15,i)); %for large freq range
end

    % 1st Frequency Band:
for i = 1:length(See1)
    [L1(i),~] = LDR(1/ff(8),Depth1(i)); %use ff(8) because it is the middle of the 3 freqs
    c1(i) = L1(i)/(1/ff(8));
    Cg1(i) = (c1(i)/2)*(1 + ((4*pi*Depth1(i)/L1(i))/sinh(4*pi*Depth1(i)/L1(i))));

    [L2(i),~] = LDR(1/ff(8),Depth2(i)); %use ff(8) because it is the middle of the 3 freqs
    c2(i) = L2(i)/(1/ff(8));
    Cg2(i) = (c2(i)/2)*(1 + ((4*pi*Depth2(i)/L2(i))/sinh(4*pi*Depth2(i)/L2(i))));
end


    % Calculate the wave height using shoaling and ref coeff:
H_SR_snell{2} = H_SR{2}./(sqrt(Cg2./Cg1).*sqrt(cosd(EWdir_SR1)./cosd(EWdir_SR2))); %snell height at B05 (SR)


% Snell:
theta2_1stBand{2} = asind((Cg2./Cg1).*sind(EWdir_SR1));
Compare_1stBand{2} = meanangle(ObsTheta2_SR);
% Plot:
figure(3);clf; 
plot(theta2_1stBand{2},'r','LineWidth',2); hold on; grid on;
plot(Compare_1stBand{2},'b','LineWidth',2);
title('B05: Snells Law vs Observed (Freq: 0.0686-0.0882 Hz)')
legend('Snells Dir','Obs Dir')
ylabel('Wave Dir (relative to shore normal)')



    % 2nd Frequency Band:
for i = 1:length(See1)
    [L1(i),~] = LDR(1/ff(13),Depth1(i)); %use ff(13) because it is the middle of the 5 freqs
    c1(i) = L1(i)/(1/ff(13));
    Cg1(i) = (c1(i)/2)*(1 + ((4*pi*Depth1(i)/L1(i))/sinh(4*pi*Depth1(i)/L1(i))));

    [L2(i),~] = LDR(1/ff(13),Depth2(i)); %use ff(13) because it is the middle of the 5 freqs
    c2(i) = L2(i)/(1/ff(13));
    Cg2(i) = (c2(i)/2)*(1 + ((4*pi*Depth2(i)/L2(i))/sinh(4*pi*Depth2(i)/L2(i))));
end


    % Calculate the wave height using shoaling and ref coeff:
H_LR_snell{2} = H_LR{2}./(sqrt(Cg2./Cg1).*sqrt(cosd(EWdir_LR1)./cosd(EWdir_LR2))); %snell height at B05 (LR)


% Snell:
theta2_2ndBand{2} = asind((Cg2./Cg1).*sind(EWdir_LR1));
Compare_2ndBand{2} = meanangle(ObsTheta2_LR);
% Plot:
figure(4);clf; 
plot(theta2_2ndBand{2},'r','LineWidth',2); hold on; grid on;
plot(Compare_2ndBand{2},'b','LineWidth',2);
title('B05: Snells Law vs Observed (Freq: 0.1078-0.1470 Hz)')
legend('Snells Dir','Obs Dir')
ylabel('Wave Dir (relative to shore normal)')


%% Plot the observed wave heights against snell wave heights

figure(5);clf;

subplot(2,2,1)
plot(H_SR{2}([1:397 505:end]))
hold on; grid on;
plot(H_SR_snell{1})
title('Wave Height: B03 SR')
ylabel('H (m)')
legend('Observed (mean=0.201m)','Snell (mean=0.195m)')

subplot(2,2,2)
plot(H_LR{2}([1:397 505:end]))
hold on; grid on;
plot(H_LR_snell{1})
title('Wave Height: B03 LR')
ylabel('H (m)')
legend('Observed (mean=0.324m)','Snell (mean=0.332m)')

subplot(2,2,3)
plot(H_SR{3})
hold on; grid on;
plot(H_SR_snell{2})
title('Wave Height: B05 SR')
ylabel('H (m)')
legend('Observed (mean=0.211m)','Snell (mean=0.208m)')

subplot(2,2,4)
plot(H_LR{3})
hold on; grid on;
plot(H_LR_snell{2})
title('Wave Height: B05 LR')
ylabel('H (m)')
legend('Observed (mean=0.295m)','Snell (mean=0.305m)')



%% Save Variables

Snell_SR = theta2_1stBand;
Snell_LR = theta2_2ndBand;

Obs_SR = Compare_1stBand;
Obs_LR = Compare_2ndBand;

save('Snells.mat','Snell_SR','Snell_LR','Obs_SR','Obs_LR','B01obs_SR')


% SR - 36 32
% LR - 11 3














%%
%     % In the 1st frequency band:
%     
% for i = 1:3
%     % Determine group speed:
%     ind = i+6; %this is based on the indminf1 above
%     [L1(i),k1(i)] = LDR(1/ff(ind),BmeanDepth1);
%     c = L1(i)/(1/ff(i));
%     Cg1(i) = (c/2)*(1 + ((4*pi*BmeanDepth1/L1(i))/sinh(4*pi*BmeanDepth1/L1(i))));
%     
%     [L2(i),k2(i)] = LDR(1/ff(ind),BmeanDepth2);
%     c = L2(i)/(1/ff(i));
%     Cg2(i) = (c/2)*(1 + ((4*pi*BmeanDepth2/L2(i))/sinh(4*pi*BmeanDepth2/L2(i))));
%     
%     %Snell:
%     theta2(i) = asind((Cg2(i)/Cg1(i))*sind(ObsTheta1(ind)));
% end
% 
% disp('1st Frequency Band (0.0686 - 0.0882 Hz):')
% theta2
% compare = ObsTheta2(indminf1:indmaxf1)'
% 
% 
% 
% % In the 2nd frequency band:
% 
% for i = 1:5
%     % Determine group speed:
%     ind = i+10; %this is based on the indminf2 above
%     [L1(i),k1(i)] = LDR(1/ff(ind),BmeanDepth1);
%     c = L1(i)/(1/ff(i));
%     Cg1(i) = (c/2)*(1 + ((4*pi*BmeanDepth1/L1(i))/sinh(4*pi*BmeanDepth1/L1(i))));
%     
%     [L2(i),k2(i)] = LDR(1/ff(ind),BmeanDepth2);
%     c = L2(i)/(1/ff(i));
%     Cg2(i) = (c/2)*(1 + ((4*pi*BmeanDepth2/L2(i))/sinh(4*pi*BmeanDepth2/L2(i))));
%     
%     %Snell:
%     theta2(i) = asind((Cg2(i)/Cg1(i))*sind(ObsTheta1(ind)));
% end
% 
% disp('2nd Frequency Band (0.1078 - 0.1470 Hz):')
% theta2
% compare = ObsTheta2(indminf2:indmaxf2)'
% 



