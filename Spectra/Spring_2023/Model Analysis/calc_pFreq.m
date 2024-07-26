%% calc_pFreq.m
%
% Noah Clark        5/16/24
%
% Purpose: - To calculate the peak frequency for each hour at each buoy for
%            the observed and model (SWAN050, SWAN050sm, and COUP050) data 
%          - Plot the model peak frequency results against the observed
%            peak frequencies (3 plots)
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

%% Preliminaries

clc;clear;

load('InterpModelObs.mat')

addpath('../Stat Functions') % the functions for the statistical analysis


%% Determine Peak Frequencies

Bnames = fields(Buoy.Obs); % to list buoy names

    % Observed:
for i = 1:9
    eval(sprintf('STODO = Buoy.Obs.%s;',Bnames{i}))

    for j = 1:length(Buoy.time)
        [~,ind_pf] = max(STODO.See(:,j));
        pf(j) = Buoy.freq(ind_pf);
    end

    eval(sprintf('Buoy.Obs.%s.pfreq = pf;',Bnames{i}))
end


    % SWAN050:
for i = 1:9
    eval(sprintf('STODO = Buoy.SWAN050.%s;',Bnames{i}))

    for j = 1:length(Buoy.time)
        [~,ind_pf] = max(STODO.See(:,j));
        pf(j) = Buoy.freq(ind_pf);
    end

    eval(sprintf('Buoy.SWAN050.%s.pfreq = pf;',Bnames{i}))
end


    % SWAN050sm:
for i = 1:9
    eval(sprintf('STODO = Buoy.SWAN050sm.%s;',Bnames{i}))

    for j = 1:length(Buoy.time)
        [~,ind_pf] = max(STODO.See(:,j));
        pf(j) = Buoy.freq(ind_pf);
    end

    eval(sprintf('Buoy.SWAN050sm.%s.pfreq = pf;',Bnames{i}))
end


    % COUP050:
for i = 1:9
    eval(sprintf('STODO = Buoy.COUP050.%s;',Bnames{i}))

    for j = 1:length(Buoy.time)
        [~,ind_pf] = max(STODO.See(:,j));
        pf(j) = Buoy.freq(ind_pf);
    end

    eval(sprintf('Buoy.COUP050.%s.pfreq = pf;',Bnames{i}))
end


    % SWAN075:
for i = 1:9
    eval(sprintf('STODO = Buoy.SWAN075.%s;',Bnames{i}))

    for j = 1:length(Buoy.time)
        [~,ind_pf] = max(STODO.See(:,j));
        pf(j) = Buoy.freq(ind_pf);
    end

    eval(sprintf('Buoy.SWAN075.%s.pfreq = pf;',Bnames{i}))
end


%% Plot

pformat = {'xb','xm','xr','xc','xk','ob','om','or','oc'};

one_to_one = 0.04:0.01:0.27;

CATO = [];
CATM{1} = []; CATM{2} = []; CATM{3} = []; CATM{4} = [];

figure(1);clf;

for i = 1:9
    
    eval(sprintf('pfobs = Buoy.Obs.%s.pfreq;',Bnames{i}))
    CATO = cat(2,CATO,pfobs);
    
    eval(sprintf('pfM = Buoy.SWAN050.%s.pfreq;',Bnames{i}))
    subplot(141)
    plot(pfobs,pfM,pformat{i},'MarkerSize',8,'LineWidth',2)
    hold on
    CATM{1} = cat(2,CATM{1},pfM);
    
    eval(sprintf('pfM = Buoy.SWAN050sm.%s.pfreq;',Bnames{i}))
    subplot(142)
    plot(pfobs,pfM,pformat{i},'MarkerSize',8,'LineWidth',2)
    hold on
    CATM{2} = cat(2,CATM{2},pfM);

    eval(sprintf('pfM = Buoy.COUP050.%s.pfreq;',Bnames{i}))
    subplot(143)
    plot(pfobs,pfM,pformat{i},'MarkerSize',8,'LineWidth',2)
    hold on
    CATM{3} = cat(2,CATM{3},pfM);
    
    eval(sprintf('pfM = Buoy.SWAN075.%s.pfreq;',Bnames{i}))
    subplot(144)
    plot(pfobs,pfM,pformat{i},'MarkerSize',8,'LineWidth',2)
    hold on
    CATM{4} = cat(2,CATM{4},pfM);

   
end

subplot(141)
plot(one_to_one,one_to_one,'--k')
axis equal; grid on;
legend(Bnames,'NumColumns',2,'fontsize',12)
xlabel('Obs f_p (Hz)','fontsize',13)
ylabel('SWAN050 f_p (Hz)','fontsize',13)
title('SWAN050 f_p Comparison','fontsize',16)
[~,Sk,maxcor,rmse] = mod_error(CATM{1},CATO);
str = {sprintf('RMSE = %.4f',round(rmse.r1,4)),...
    sprintf('Sk = %.4f',round(Sk.sk1,4)),...
    sprintf('CC = %.4f',round(maxcor.CCall(101),4))};
dim = [0.161405497168987,0.146245096729904,0.132255386023676,0.152416352888909];
annotation('textbox',dim,'String',str,'fontsize',14,'EdgeColor','none')

subplot(142)
plot(one_to_one,one_to_one,'--k')
axis equal; grid on;
xlabel('Obs f_p (Hz)','fontsize',13)
ylabel('SWAN050sm f_p (Hz)','fontsize',13)
title('SWAN050sm f_p Comparison','fontsize',16)
[~,Sk,maxcor,rmse] = mod_error(CATM{2},CATO);
str = {sprintf('RMSE = %.4f',round(rmse.r1,4)),...
    sprintf('Sk = %.4f',round(Sk.sk1,4)),...
    sprintf('CC = %.4f',round(maxcor.CCall(101),4))};
dim = [0.36665832572836,0.142205356887847,0.132255386023676,0.152416352888909];
annotation('textbox',dim,'String',str,'fontsize',14,'EdgeColor','none')

subplot(143)
plot(one_to_one,one_to_one,'--k')
axis equal; grid on;
xlabel('Obs f_p (Hz)','fontsize',13)
ylabel('COUP050 f_p (Hz)','fontsize',13)
title('COUP050 f_p Comparison','fontsize',16)
[~,Sk,maxcor,rmse] = mod_error(CATM{3},CATO);
str = {sprintf('RMSE = %.4f',round(rmse.r1,4)),...
    sprintf('Sk = %.4f',round(Sk.sk1,4)),...
    sprintf('CC = %.4f',round(maxcor.CCall(101),4))};
dim = [0.573546223517521,0.151069346887781,0.132255386023676,0.152416352888909];
annotation('textbox',dim,'String',str,'fontsize',14,'EdgeColor','none')

subplot(144)
plot(one_to_one,one_to_one,'--k')
axis equal; grid on;
xlabel('Obs f_p (Hz)','fontsize',13)
ylabel('SWAN075 f_p (Hz)','fontsize',13)
title('SWAN075 f_p Comparison','fontsize',16)
[~,Sk,maxcor,rmse] = mod_error(CATM{4},CATO);
str = {sprintf('RMSE = %.4f',round(rmse.r1,4)),...
    sprintf('Sk = %.4f',round(Sk.sk1,4)),...
    sprintf('CC = %.4f',round(maxcor.CCall(101),4))};
dim = [0.781386460203912,0.142312779462212,0.132255386023676,0.152416352888909];
annotation('textbox',dim,'String',str,'fontsize',14,'EdgeColor','none')


%%

clear i j ind_pf pf STODO

save('InterpModelObs.mat','-append','Buoy')






