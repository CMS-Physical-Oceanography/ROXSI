%% calc_ShoalCoeff.m
%
% Noah Clark            1/27/23
%
%
% Purpose: To calculate the shoaling coefficients between the instruments 
%          by E1/E2 and compare it against the shoaling coefficient using
%          sqrt(Cg2/Cg1) at the three time periods
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

%% Preliminaries

clc;clear;

load('3Periods.mat')
load('WBvariables.mat','Bdepth')

ff = 0.0098:0.0098:(0.0098*129);
ffs = 0.0098:0.0098:(0.0098*35); %a shortened version of ff (for interpolated ADCP data)
ffa = 0.0028:0.0028:(0.0028*126); %freq vector for original ADCP data
TT = 1./ff;
TTs = 1./ffs;

Def_Freq = 0.088; % >: Sea      <: Swell

    % ADCP Data:
addpath('../Fall (0825)/ADCP Data')
% Load data for B10
load('roxsi_signature_L2_B10_See.mat')
See_B10 = A.See;
See_B10 = interp1(ffa,See_B10,ff);
See_B10 = See_B10(1:35,:);
SeaP.B10.mSee = mean(See_B10(:,SeaP.ind1-151:SeaP.ind2-151),2);
SwellP.B10.mSee = mean(See_B10(:,SwellP.ind1-151:SwellP.ind2-151),2);
BothP.B10.mSee = mean(See_B10(:,BothP.ind1-151:BothP.ind2-151),2);
Depth_B10 = A.depth;
% Load data for B13
load('roxsi_signature_L2_B13_See.mat')
See_B13 = A.See;
See_B13 = interp1(ffa,See_B13,ff);
See_B13 = See_B13(1:35,:);
SeaP.B13.mSee = mean(See_B13(:,SeaP.ind1-151:SeaP.ind2-151),2);
SwellP.B13.mSee = mean(See_B13(:,SwellP.ind1-151:SwellP.ind2-151),2);
BothP.B13.mSee = mean(See_B13(:,BothP.ind1-151:BothP.ind2-151),2);
Depth_B13 = A.depth;

clear See_B10 See_B13 


%% Sea Period: (B01-B03, B03-B05, B05-B10, B10-B13)

    % Calculate:
i1 = SeaP.ind1;
i2 = SeaP.ind2;
% Depths
SeaP.B01.avgD = mean(Bdepth{1}(i1-504:i2-504));
SeaP.B03.avgD = mean(Bdepth{2}(i1:i2));
SeaP.B05.avgD = mean(Bdepth{3}(i1:i2));
SeaP.B10.avgD = mean(Depth_B10(i1-151:i2-151));
SeaP.B13.avgD = mean(Depth_B13(i1-151:i2-151));

% Wavelength
for i = 1:129
    [SeaP.B01.L(i),~] = LDR(TT(i),SeaP.B01.avgD);
    [SeaP.B03.L(i),~] = LDR(TT(i),SeaP.B03.avgD);
    [SeaP.B05.L(i),~] = LDR(TT(i),SeaP.B05.avgD);
end
for i = 1:35
    [SeaP.B10.L(i),~] = LDR(TTs(i),SeaP.B10.avgD);
    [SeaP.B13.L(i),~] = LDR(TTs(i),SeaP.B13.avgD);
end
% Speed
SeaP.B01.c = SeaP.B01.L./TT;
SeaP.B03.c = SeaP.B03.L./TT;
SeaP.B05.c = SeaP.B05.L./TT;
SeaP.B10.c = SeaP.B10.L./TTs;
SeaP.B13.c = SeaP.B13.L./TTs;

% Group Speed
SeaP.B01.Cg = (SeaP.B01.c./2).*(1+((4*pi*SeaP.B01.avgD)./SeaP.B01.L)./sinh((4*pi*SeaP.B01.avgD)./SeaP.B01.L));
SeaP.B03.Cg = (SeaP.B03.c./2).*(1+((4*pi*SeaP.B03.avgD)./SeaP.B03.L)./sinh((4*pi*SeaP.B03.avgD)./SeaP.B03.L));
SeaP.B05.Cg = (SeaP.B05.c./2).*(1+((4*pi*SeaP.B05.avgD)./SeaP.B05.L)./sinh((4*pi*SeaP.B05.avgD)./SeaP.B05.L));
SeaP.B10.Cg = (SeaP.B10.c./2).*(1+((4*pi*SeaP.B10.avgD)./SeaP.B10.L)./sinh((4*pi*SeaP.B10.avgD)./SeaP.B10.L));
SeaP.B13.Cg = (SeaP.B13.c./2).*(1+((4*pi*SeaP.B13.avgD)./SeaP.B13.L)./sinh((4*pi*SeaP.B13.avgD)./SeaP.B13.L));

% Cg Shoaling Coefficient (Ks_Cg)
SeaP.Ks_Cg{1} = sqrt(SeaP.B01.Cg./SeaP.B03.Cg);
SeaP.Ks_Cg{2} = sqrt(SeaP.B03.Cg./SeaP.B05.Cg);
SeaP.Ks_Cg{3} = sqrt(SeaP.B05.Cg(1:35)./SeaP.B10.Cg);
SeaP.Ks_Cg{4} = sqrt(SeaP.B10.Cg./SeaP.B13.Cg);

% E Shoaling Coefficeint (Ks_E)
SeaP.Ks_E{1} = SeaP.B03.mSee./SeaP.B01.mSee;
SeaP.Ks_E{2} = SeaP.B05.mSee./SeaP.B03.mSee;
SeaP.Ks_E{3} = SeaP.B10.mSee./SeaP.B05.mSee(1:35);
SeaP.Ks_E{4} = SeaP.B13.mSee./SeaP.B10.mSee;


    % Plot:
figure(1);clf;

subplot(4,1,1)
plot(ff,SeaP.Ks_Cg{1},'b','LineWidth',1.5)
hold on; grid on;
plot(ff,SeaP.Ks_E{1},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Sea Period: Shoaling Coefficient from B01 to B03','fontsize',18)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0.7 1.5])

subplot(4,1,2)
plot(ff,SeaP.Ks_Cg{2},'b','LineWidth',1.5)
hold on; grid on;
plot(ff,SeaP.Ks_E{2},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Sea Period: Shoaling Coefficient from B03 to B05','fontsize',18)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0.7 1.5])

subplot(4,1,3)
plot(ffs,SeaP.Ks_Cg{3},'b','LineWidth',1.5)
hold on; grid on;
plot(ffs,SeaP.Ks_E{3},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Sea Period: Shoaling Coefficient from B05 to B10','fontsize',18)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0 9.5])

subplot(4,1,4)
plot(ffs,SeaP.Ks_Cg{4},'b','LineWidth',1.5)
hold on; grid on;
plot(ffs,SeaP.Ks_E{4},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Sea Period: Shoaling Coefficient from B10 to B13','fontsize',18)
xlabel('f [Hz]','fontsize',17)
ylabel('K_s')
xlim([0 0.34])
ylim([0 2])



%% Swell Period: (B01-B03, B03-B05, B05-B10, B10-B13)

    % Calculate
i1 = SwellP.ind1;
i2 = SwellP.ind2;
% Depths
SwellP.B01.avgD = mean(Bdepth{1}(i1:i2));
SwellP.B03.avgD = mean(Bdepth{2}(i1:i2));
SwellP.B05.avgD = mean(Bdepth{3}(i1:i2));
SwellP.B10.avgD = mean(Depth_B10(i1-151:i2-151));
SwellP.B13.avgD = mean(Depth_B13(i1-151:i2-151));

% Wavelength
for i = 1:129
    [SwellP.B01.L(i),~] = LDR(TT(i),SwellP.B01.avgD);
    [SwellP.B03.L(i),~] = LDR(TT(i),SwellP.B03.avgD);
    [SwellP.B05.L(i),~] = LDR(TT(i),SwellP.B05.avgD);
end
for i = 1:35
    [SwellP.B10.L(i),~] = LDR(TTs(i),SwellP.B10.avgD);
    [SwellP.B13.L(i),~] = LDR(TTs(i),SwellP.B13.avgD);
end

% Speed
SwellP.B01.c = SwellP.B01.L./TT;
SwellP.B03.c = SwellP.B03.L./TT;
SwellP.B05.c = SwellP.B05.L./TT;
SwellP.B10.c = SwellP.B10.L./TTs;
SwellP.B13.c = SwellP.B13.L./TTs;

% Group Speed
SwellP.B01.Cg = (SwellP.B01.c./2).*(1+((4*pi*SwellP.B01.avgD)./SwellP.B01.L)./sinh((4*pi*SwellP.B01.avgD)./SwellP.B01.L));
SwellP.B03.Cg = (SwellP.B03.c./2).*(1+((4*pi*SwellP.B03.avgD)./SwellP.B03.L)./sinh((4*pi*SwellP.B03.avgD)./SwellP.B03.L));
SwellP.B05.Cg = (SwellP.B05.c./2).*(1+((4*pi*SwellP.B05.avgD)./SwellP.B05.L)./sinh((4*pi*SwellP.B05.avgD)./SwellP.B05.L));
SwellP.B10.Cg = (SwellP.B10.c./2).*(1+((4*pi*SwellP.B10.avgD)./SwellP.B10.L)./sinh((4*pi*SwellP.B10.avgD)./SwellP.B10.L));
SwellP.B13.Cg = (SwellP.B13.c./2).*(1+((4*pi*SwellP.B13.avgD)./SwellP.B13.L)./sinh((4*pi*SwellP.B13.avgD)./SwellP.B13.L));

% Cg Shoaling Coefficient (Ks_Cg)
SwellP.Ks_Cg{1} = sqrt(SwellP.B01.Cg./SwellP.B03.Cg);
SwellP.Ks_Cg{2} = sqrt(SwellP.B03.Cg./SwellP.B05.Cg);
SwellP.Ks_Cg{3} = sqrt(SwellP.B05.Cg(1:35)./SwellP.B10.Cg);
SwellP.Ks_Cg{4} = sqrt(SwellP.B10.Cg./SwellP.B13.Cg);

% E Shoaling Coefficeint (Ks_E)
SwellP.Ks_E{1} = SwellP.B03.mSee./SwellP.B01.mSee;
SwellP.Ks_E{2} = SwellP.B05.mSee./SwellP.B03.mSee;
SwellP.Ks_E{3} = SwellP.B10.mSee./SwellP.B05.mSee(1:35);
SwellP.Ks_E{4} = SwellP.B13.mSee./SwellP.B10.mSee;


    % Plot:
figure(2);clf;

subplot(4,1,1)
plot(ff,SwellP.Ks_Cg{1},'b','LineWidth',1.5)
hold on; grid on;
plot(ff,SwellP.Ks_E{1},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Swell Period: Shoaling Coefficient from B01 to B03','fontsize',18)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0.8 1.5])

subplot(4,1,2)
plot(ff,SwellP.Ks_Cg{2},'b','LineWidth',1.5)
hold on; grid on;
plot(ff,SwellP.Ks_E{2},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Swell Period: Shoaling Coefficient from B03 to B05','fontsize',18)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0.8 1.5])

subplot(4,1,3)
plot(ffs,SwellP.Ks_Cg{3},'b','LineWidth',1.5)
hold on; grid on;
plot(ffs,SwellP.Ks_E{3},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Swell Period: Shoaling Coefficient from B05 to B10','fontsize',18)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0 3])

subplot(4,1,4)
plot(ffs,SwellP.Ks_Cg{4},'b','LineWidth',1.5)
hold on; grid on;
plot(ffs,SwellP.Ks_E{4},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Swell Period: Shoaling Coefficient from B10 to B13','fontsize',18)
xlabel('f [Hz]','fontsize',17)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0 3])


%% Both Period: (B01-B03, B03-B05, B05-B10, B10-B13)

    % Calculate:
i1 = BothP.ind1;
i2 = BothP.ind2;
% Depths
BothP.B01.avgD = mean(Bdepth{1}(i1:i2));
BothP.B03.avgD = mean(Bdepth{2}(i1:i2));
BothP.B05.avgD = mean(Bdepth{3}(i1:i2));
BothP.B10.avgD = mean(Depth_B10(i1-151:i2-151));
BothP.B13.avgD = mean(Depth_B13(i1-151:i2-151));

% Wavelength
for i = 1:129
    [BothP.B01.L(i),~] = LDR(TT(i),BothP.B01.avgD);
    [BothP.B03.L(i),~] = LDR(TT(i),BothP.B03.avgD);
    [BothP.B05.L(i),~] = LDR(TT(i),BothP.B05.avgD);
end
for i = 1:35
    [BothP.B10.L(i),~] = LDR(TTs(i),BothP.B10.avgD);
    [BothP.B13.L(i),~] = LDR(TTs(i),BothP.B13.avgD);
end

% Speed
BothP.B01.c = BothP.B01.L./TT;
BothP.B03.c = BothP.B03.L./TT;
BothP.B05.c = BothP.B05.L./TT;
BothP.B10.c = BothP.B10.L./TTs;
BothP.B13.c = BothP.B13.L./TTs;

% Group Speed
BothP.B01.Cg = (BothP.B01.c./2).*(1+((4*pi*BothP.B01.avgD)./BothP.B01.L)./sinh((4*pi*BothP.B01.avgD)./BothP.B01.L));
BothP.B03.Cg = (BothP.B03.c./2).*(1+((4*pi*BothP.B03.avgD)./BothP.B03.L)./sinh((4*pi*BothP.B03.avgD)./BothP.B03.L));
BothP.B05.Cg = (BothP.B05.c./2).*(1+((4*pi*BothP.B05.avgD)./BothP.B05.L)./sinh((4*pi*BothP.B05.avgD)./BothP.B05.L));
BothP.B10.Cg = (BothP.B10.c./2).*(1+((4*pi*BothP.B10.avgD)./BothP.B10.L)./sinh((4*pi*BothP.B10.avgD)./BothP.B10.L));
BothP.B13.Cg = (BothP.B13.c./2).*(1+((4*pi*BothP.B13.avgD)./BothP.B13.L)./sinh((4*pi*BothP.B13.avgD)./BothP.B13.L));

% Cg Shoaling Coefficient (Ks_Cg)
BothP.Ks_Cg{1} = sqrt(BothP.B01.Cg./BothP.B03.Cg);
BothP.Ks_Cg{2} = sqrt(BothP.B03.Cg./BothP.B05.Cg);
BothP.Ks_Cg{3} = sqrt(BothP.B05.Cg(1:35)./BothP.B10.Cg);
BothP.Ks_Cg{4} = sqrt(BothP.B10.Cg./BothP.B13.Cg);

% E Shoaling Coefficeint (Ks_E)
BothP.Ks_E{1} = BothP.B03.mSee./BothP.B01.mSee;
BothP.Ks_E{2} = BothP.B05.mSee./BothP.B03.mSee;
BothP.Ks_E{3} = BothP.B10.mSee./BothP.B05.mSee(1:35);
BothP.Ks_E{4} = BothP.B13.mSee./BothP.B10.mSee;


    % Plot:
figure(3);clf;

subplot(4,1,1)
plot(ff,BothP.Ks_Cg{1},'b','LineWidth',1.5)
hold on; grid on;
plot(ff,BothP.Ks_E{1},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Both Period: Shoaling Coefficient from B01 to B03','fontsize',18)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0.8 1.65])

subplot(4,1,2)
plot(ff,BothP.Ks_Cg{2},'b','LineWidth',1.5)
hold on; grid on;
plot(ff,BothP.Ks_E{2},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Both Period: Shoaling Coefficient from B03 to B05','fontsize',18)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0.8 1.65])

subplot(4,1,3)
plot(ffs,BothP.Ks_Cg{3},'b','LineWidth',1.5)
hold on; grid on;
plot(ffs,BothP.Ks_E{3},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Both Period: Shoaling Coefficient from B05 to B10','fontsize',18)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0 6.8])

subplot(4,1,4)
plot(ffs,BothP.Ks_Cg{4},'b','LineWidth',1.5)
hold on; grid on;
plot(ffs,BothP.Ks_E{4},'r','LineWidth',1.5)
xline(Def_Freq,'--m')
legend('K_s: C_g','K_s: E','Sea vs Swell','fontsize',12)
title('Both Period: Shoaling Coefficient from B10 to B13','fontsize',18)
xlabel('f [Hz]','fontsize',17)
ylabel('K_s','fontsize',16)
xlim([0 0.34])
ylim([0 2])


%% Create Matricies of Ks_Cg
% - The Ks_Cg will be the same for Sea, Swell, and Both
% - Just use the SeaP Ks_Cgs

Matrix_Ks_Cg = [SeaP.Ks_Cg{1}(1:21);SeaP.Ks_Cg{2}(1:21);...
    SeaP.Ks_Cg{3}(1:21);SeaP.Ks_Cg{4}(1:21)];


%% Figure for OSM
% Find same plot but for entire time average 

% Average Sees
load('WBvariables.mat','BSee')
avgSee1 = mean(BSee{1},2);
avgSee2 = mean(BSee{2}(:,[1:397 505:end]),2);
% Depth
avgD1 = mean(Bdepth{1});
avgD2 = mean(Bdepth{2}([1:397 505:end]));
% Wavelength
for i = 1:129
    [L1(i),~] = LDR(TT(i),avgD1);
    [L2(i),~] = LDR(TT(i),avgD2);
end
% Speed
c1 = L1./TT;
c2 = L2./TT;
% Group Speed
Cg1 = (c1./2).*(1+((4*pi*avgD1)./L1)./sinh((4*pi*avgD1)./L1));
Cg2 = (c2./2).*(1+((4*pi*avgD2)./L2)./sinh((4*pi*avgD2)./L2));
% Cg Shoaling Coefficient (Ks_Cg)
Ks_Cg = sqrt(Cg1./Cg2);
% E Shoaling Coefficeint (Ks_E)
Ks_E = avgSee2./avgSee1;


figure(4);clf;
plot(ff,Ks_Cg,'b','LineWidth',2)
hold on; grid on;
plot(ff,Ks_E,'r','LineWidth',2)
xline(Def_Freq,'--k') % line for sea vs swell
yline(1,'k','LineWidth',1.2) % line at 1 to differentiate <1 and <1 Ks
ax = gca;
ax.FontSize = 15;
legend('Theoretical','Observed','Sea vs Swell','fontsize',19)
title('Shoaling Coefficient from B01 to B03','fontsize',22)
ylabel('K_s','fontsize',19)
xlabel('Frequency (Hz)','fontsize',19)
xlim([0 0.35])
ylim([0.8 2])


%% Save Variables and Figures

clear i i1 i2 A

% Save Variables
save('3Periods.mat','SeaP','SwellP','BothP')

% Save Figures
cd 'Spring Figures/Shoal_Coeff_Plots'
saveas(figure(1),'Sea_ShoalCoef.jpeg')
saveas(figure(2),'Swell_ShoalCoef.jpeg')
saveas(figure(3),'Both_ShoalCoef.jpeg')






