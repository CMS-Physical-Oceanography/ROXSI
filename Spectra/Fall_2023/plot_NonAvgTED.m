%% plot_NonAvgTED.m
%
% Noah Clark
% Date Created: 9/29/23
%
% Purpose: 
%         Determine and plot the Nielsen and observed dissipation using the
%         functions NC_NeilsonDiss and NC_ObsDiss at each time and
%         frequency. This is done between all buoys and ADCPs to determine
%         how similar they are at all points in time and frequencies by
%         plotting spectrograms of these dissipations.
%
%
% TO-DO:
%       - can still play around with caxis to make figures look better
%
%

%% Preliminaries

clc;clear;
load('SM&ADCP_All.mat')

%addpath('C:\Users\nsc4750\OneDrive - UNC-Wilmington\CMS Work\Summer (0516-0728)')
load('WBvariables.mat','Bmeanspec','Xmeanspec','Butm','Xutm','BSee',...
    'XSee','BEMEM','XEMEM','Xdepth','Bdepth','Btime','Xtime',...
    'BNormWaveDir','XNormWaveDir')

%FunctionPath = dir('/Users/noahclark/Library/CloudStorage/OneDrive-UNC-Wilmington/CMS Work/Summer (0516-0728)/Main Functions');
%addpath(FunctionPath.folder)

addpath('../Summer (0516-0728)/Main Functions')

a1 = 5.5; a2 = -0.2; a3 = -6.3;
kw = 1;
ff = (0.0098:0.0098:0.343);
df = 0.0098;


%% X01 to X03

See1 = XSee{1}(1:35,[1:172 174:end]);
See2 = XSee{2}(1:35,[1:172 174:end]);
Direc1 = XEMEM{1}(1:35,[1:172 174:end]);
Direc2 = XEMEM{2}(1:35,[1:172 174:end]);
utm1 = Xutm{1};
utm2 = Xutm{2};
Depth1 = Xdepth{1}([1:172 174:end]);
Depth2 = Xdepth{2}([1:172 174:end]);
loop = size(See1);
Time = Xtime{1}([1:172 174:end]);

for j = 1:loop(2)
    Obs_TED{1}(:,j) = NC_ObsDiss(See1(:,j),See2(:,j),Direc1(:,j),...
        Direc2(:,j),ff,utm1,utm2,Depth1(j),Depth2(j));
    [~,~,~,N_TED1(:,j)] = NC_NeilsonDiss(See1(:,j),ff,Depth1(j),...
        kw,a1,a2,a3,'r');
    [~,~,~,N_TED2(:,j)] = NC_NeilsonDiss(See2(:,j),ff,Depth2(j),...
        kw,a1,a2,a3,'r');
end

N_TED{1} = (N_TED1 + N_TED2)./2;



figure(1);clf;
subplot(2,1,1)
pcolor(Time,ff,Obs_TED{1})
cb = colorbar;
shading interp
caxis([-2 5])
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Observed Dissipation from X01 to X03')

subplot(2,1,2)
pcolor(Time,ff,N_TED{1})
cb = colorbar;
shading interp
caxis([-2 5])
xlabel('Date')
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Nielson Dissipation from X01 to X03')




%% X03 to X04

See1 = XSee{2}(1:35,[1:172 174:end]);
See2 = XSee{3}(1:35,[1:172 174:end]);
Direc1 = XEMEM{2}(1:35,[1:172 174:end]);
Direc2 = XEMEM{3}(1:35,[1:172 174:end]);
utm1 = Xutm{2};
utm2 = Xutm{3};
Depth1 = Xdepth{2}([1:172 174:end]);
Depth2 = Xdepth{3}([1:172 174:end]);
loop = size(See2);
Time = Xtime{2}([1:172 174:end]);

clear N_TED1 N_TED2

for j = 1:loop(2)
    Obs_TED{2}(:,j) = NC_ObsDiss(See1(:,j),See2(:,j),Direc1(:,j),...
        Direc2(:,j),ff,utm1,utm2,Depth1(j),Depth2(j));
    [~,~,~,N_TED1(:,j)] = NC_NeilsonDiss(See1(:,j),ff,Depth1(j),...
        kw,a1,a2,a3,'r');
    [~,~,~,N_TED2(:,j)] = NC_NeilsonDiss(See2(:,j),ff,Depth2(j),...
        kw,a1,a2,a3,'r');
end
N_TED{2} = (N_TED1 + N_TED2)./2;

figure(2);clf;
subplot(2,1,1)
pcolor(Time,ff,Obs_TED{2})
cb = colorbar;
shading interp
caxis([-2 5])
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Observed Dissipation from X03 to X04')

subplot(2,1,2)
pcolor(Time,ff,N_TED{2})
cb = colorbar;
shading interp
caxis([-2 5])
xlabel('Date')
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Nielson Dissipation from X03 to X04')


%% X04 to X05

See1 = XSee{3}(1:35,152:end);
See2 = interp1(ADCP.freq,ADCP.X05.See(:,1:681),ff);

Direc1 = XEMEM{3}(1:35,152:end);
Direc2 = Direc1;
utm1 = Xutm{3};
utm2 = ADCP.X05.utm;
Depth1 = Xdepth{3}(152:end);
Depth2 = ADCP.X05.depth(1:681);
loop = size(See2);
Time = Xtime{3}(152:end);

clear N_TED1 N_TED2

for j = 1:loop(2)
    Obs_TED{3}(:,j) = NC_ObsDiss(See1(:,j),See2(:,j),Direc1(:,j),...
        Direc2(:,j),ff,utm1,utm2,Depth1(j),Depth2(j));
    [~,~,~,N_TED1(:,j)] = NC_NeilsonDiss(See1(:,j),ff,Depth1(j),...
        kw,a1,a2,a3,'r');
    [~,~,~,N_TED2(:,j)] = NC_NeilsonDiss(See2(:,j),ff,Depth2(j),...
        kw,a1,a2,a3,'r');
end
N_TED{3} = (N_TED1 + N_TED2)./2;

figure(3);clf;
subplot(2,1,1)
pcolor(Time,ff,Obs_TED{3})
cb = colorbar;
shading interp
caxis([-2 5])
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Observed Dissipation from X04 to X05')

subplot(2,1,2)
pcolor(Time,ff,N_TED{3})
cb = colorbar;
shading interp
caxis([-2 5])
xlabel('Date')
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Nielson Dissipation from X04 to X05')


%% From X05 to X11

See1 = interp1(ADCP.freq,ADCP.X05.See(:,1:681),ff);
See2 = interp1(ADCP.freq,ADCP.X11.See(:,1:681),ff);
Direc1 = XNormWaveDir.*ones(35,681);
Direc2 = Direc1;
utm1 = ADCP.X05.utm;
utm2 = ADCP.X11.utm;
Depth1 = ADCP.X05.depth(1:681);
Depth2 = ADCP.X11.depth(1:681);
loop = size(See1);
Time = Xtime{3}(152:end);

clear N_TED1 N_TED2

for j = 1:loop(2)
    Obs_TED{4}(:,j) = NC_ObsDiss(See1(:,j),See2(:,j),Direc1(:,j),...
        Direc2(:,j),ff,utm1,utm2,Depth1(j),Depth2(j));
    [~,~,~,N_TED1(:,j)] = NC_NeilsonDiss(See1(:,j),ff,Depth1(j),...
        kw,a1,a2,a3,'r');
    [~,~,~,N_TED2(:,j)] = NC_NeilsonDiss(See2(:,j),ff,Depth2(j),...
        kw,a1,a2,a3,'r');
end

N_TED{4} = (N_TED1 + N_TED2)./2;

figure(4);clf;
subplot(2,1,1)
pcolor(Time,ff,Obs_TED{4})
cb = colorbar;
shading interp
caxis([-2 5])
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Observed Dissipation from X05 to X11')

subplot(2,1,2)
pcolor(Time,ff,N_TED{4})
cb = colorbar;
shading interp
caxis([-2 5])
xlabel('Date')
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Nielson Dissipation from X05 to X11')


%% From B01 to B03

See1 = BSee{1}(1:35,:);
See2 = BSee{2}(1:35,[1:397 505:end]);
Direc1 = BEMEM{1}(1:35,:);
Direc2 = BEMEM{2}(1:35,[1:397 505:end]);
utm1 = Butm{1};
utm2 = Butm{2};
Depth1 = Bdepth{1}(:);
Depth2 = Bdepth{2}([1:397 505:end]);
loop = size(See1);
Time = Btime{1}(:);

for j = 1:loop(2)
    Obs_TED{5}(:,j) = NC_ObsDiss(See1(:,j),See2(:,j),Direc1(:,j),...
        Direc2(:,j),ff,utm1,utm2,Depth1(j),Depth2(j));
    [~,~,~,N_TED1(:,j)] = NC_NeilsonDiss(See1(:,j),ff,Depth1(j),...
        kw,a1,a2,a3,'r');
    [~,~,~,N_TED2(:,j)] = NC_NeilsonDiss(See2(:,j),ff,Depth2(j),...
        kw,a1,a2,a3,'r');
end

N_TED{5} = (N_TED1 + N_TED2)./2;

figure(5);clf;
subplot(2,1,1)
pcolor(Time,ff,Obs_TED{5})
cb = colorbar;
shading interp
caxis([-2 5])
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Observed Dissipation from B01 to B03')

subplot(2,1,2)
pcolor(Time,ff,N_TED{5})
cb = colorbar;
shading interp
caxis([-2 5])
xlabel('Date')
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Nielson Dissipation from B01 to B03')


%% From B03 to B05

See1 = BSee{2}(1:35,:);
See2 = BSee{3}(1:35,:);
Direc1 = BEMEM{2}(1:35,:);
Direc2 = BEMEM{3}(1:35,:);
utm1 = Butm{2};
utm2 = Butm{3};
Depth1 = Bdepth{2}(:);
Depth2 = Bdepth{3}(:);
loop = size(See1);
Time = Btime{2}(:);

for j = 1:loop(2)
    Obs_TED{6}(:,j) = NC_ObsDiss(See1(:,j),See2(:,j),Direc1(:,j),...
        Direc2(:,j),ff,utm1,utm2,Depth1(j),Depth2(j));
    [~,~,~,N_TED1(:,j)] = NC_NeilsonDiss(See1(:,j),ff,Depth1(j),...
        kw,a1,a2,a3,'r');
    [~,~,~,N_TED2(:,j)] = NC_NeilsonDiss(See2(:,j),ff,Depth2(j),...
        kw,a1,a2,a3,'r');
end

N_TED{6} = (N_TED1 + N_TED2)./2;

figure(6);clf;
subplot(2,1,1)
pcolor(Time,ff,Obs_TED{6})
cb = colorbar;
shading interp
caxis([-2 5])
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Observed Dissipation from B03 to B05')

subplot(2,1,2)
pcolor(Time,ff,N_TED{6})
cb = colorbar;
shading interp
caxis([-2 5])
xlabel('Date')
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Nielson Dissipation from B03 to B05')


%% From B05 to B10

See1 = BSee{3}(1:35,152:end);
See2 = interp1(ADCP.freq,ADCP.B10.See(:,1:681),ff);
Direc1 = BEMEM{3}(1:35,152:end);
Direc2 = Direc1;
utm1 = Butm{3};
utm2 = ADCP.B10.utm;
Depth1 = Xdepth{3}(152:end);
Depth2 = ADCP.B10.depth(1:681);
loop = size(See1);
Time = Btime{3}(152:end);

clear N_TED1 N_TED2

for j = 1:loop(2)
    Obs_TED{7}(:,j) = NC_ObsDiss(See1(:,j),See2(:,j),Direc1(:,j),...
        Direc2(:,j),ff,utm1,utm2,Depth1(j),Depth2(j));
    [~,~,~,N_TED1(:,j)] = NC_NeilsonDiss(See1(:,j),ff,Depth1(j),...
        kw,a1,a2,a3,'r');
    [~,~,~,N_TED2(:,j)] = NC_NeilsonDiss(See2(:,j),ff,Depth2(j),...
        kw,a1,a2,a3,'r');
end

N_TED{7} = (N_TED1 + N_TED2)./2;

figure(7);clf;
subplot(2,1,1)
pcolor(Time,ff,Obs_TED{7})
cb = colorbar;
shading interp
caxis([-2 5])
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Observed Dissipation from B05 to B10')

subplot(2,1,2)
pcolor(Time,ff,N_TED{7})
cb = colorbar;
shading interp
caxis([-2 5])
xlabel('Date')
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Nielson Dissipation from B05 to B10')


%% From B10 to B13

See1 = interp1(ADCP.freq,ADCP.B10.See(:,1:681),ff);
See2 = interp1(ADCP.freq,ADCP.B13.See(:,1:681),ff);
Direc1 = BNormWaveDir.*ones(35,681);
Direc2 = Direc1;
utm1 = ADCP.B10.utm;
utm2 = ADCP.B13.utm;
Depth1 = ADCP.B10.depth(1:681);
Depth2 = ADCP.B13.depth(1:681);
loop = size(See1);
Time = Btime{3}(152:end);

clear N_TED1 N_TED2

for j = 1:loop(2)
    Obs_TED{8}(:,j) = NC_ObsDiss(See1(:,j),See2(:,j),Direc1(:,j),...
        Direc2(:,j),ff,utm1,utm2,Depth1(j),Depth2(j));
    [~,~,~,N_TED1(:,j)] = NC_NeilsonDiss(See1(:,j),ff,Depth1(j),...
        kw,a1,a2,a3,'r');
    [~,~,~,N_TED2(:,j)] = NC_NeilsonDiss(See2(:,j),ff,Depth2(j),...
        kw,a1,a2,a3,'r');
end

N_TED{8} = (N_TED1 + N_TED2)./2;

figure(8);clf;
subplot(2,1,1)
pcolor(Time,ff,Obs_TED{8})
cb = colorbar;
shading interp
caxis([-2 5])
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Observed Dissipation from B10 to B13')

subplot(2,1,2)
pcolor(Time,ff,N_TED{8})
cb = colorbar;
shading interp
caxis([-2 5])
xlabel('Date')
ylabel('Freq (Hz)')
ylabel(cb,'Energy Dissipation(W/m)')
title('Nielson Dissipation from B10 to B13')


%% Save Variables

obsTED = Obs_TED;
NielsonTED = N_TED;

%save('TED.mat','-append','obsTED','NielsonTED')




