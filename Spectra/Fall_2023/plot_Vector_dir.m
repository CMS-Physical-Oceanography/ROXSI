%% plot_Vector_dir.m
%
% Noah Clark        11/25/23
%
%
% Purpose: To plot the positions of the China Rock buoys in utm coordinates
%           and then plot vectors of the true wave direction in one color
%           and then the direction that is accounted for by Snell's Law in
%           another color. These Snell's directions are calculated in the
%           script calc_Snells.m.
%

clc;clear

%% Load in data

load('Snells.mat')
addpath('../Summer (0516-0728)/')
load('WBvariables.mat','Butm','BEMEM','BNormWaveDir','ShoreP')

%% Converting Direction
% Time-average the directions and convert from referencing to shore normal
%  to just regular direction (pointing in actual direction - not 180 deg
%  from)

B01obs_SR = BNormWaveDir - meanangle(B01obs_SR) - 180 - 90;

for i = 1:2
    Obs_SR{i} = BNormWaveDir - meanangle(Obs_SR{i}) - 180 - 90;
    Obs_LR{i} = BNormWaveDir - meanangle(Obs_LR{i}) - 180 - 90;
    Snell_SR{i} = BNormWaveDir - meanangle(Snell_SR{i}) - 180 - 90;
    Snell_LR{i} = BNormWaveDir - meanangle(Snell_LR{i}) - 180 - 90;
end


%% Convert angles to cartesion coordinates

B01obs_SRy = -tand(B01obs_SR);

x = 1;

for i = 1:2
    Obs_SRy{i} = -tand(Obs_SR{i});
    Obs_LRy{i} = -tand(Obs_LR{i});
    Snell_SRy{i} = -tand(Snell_SR{i});
    Snell_LRy{i} = -tand(Snell_LR{i});
end


%% Plotting

mag = 160; % the magnitude to increase the vector size


F = figure(1);clf;

quiver(Butm{1}(1),Butm{1}(2),mag*x,mag*B01obs_SRy,'r','LineWidth',2,'MarkerSize',10)
hold on; grid on;
quiver(Butm{2}(1),Butm{2}(2),mag*x,mag*Obs_SRy{1},'r','LineWidth',2,'MarkerSize',10)
quiver(Butm{2}(1),Butm{2}(2),mag*x,mag*Snell_SRy{1},'k','LineWidth',2,'MarkerSize',10)
quiver(Butm{3}(1),Butm{3}(2),mag*x,mag*Obs_SRy{2},'r','LineWidth',2,'MarkerSize',10)
quiver(Butm{3}(1),Butm{3}(2),mag*x,mag*Snell_SRy{2},'k','LineWidth',2,'MarkerSize',10)

plot(Butm{1}(1),Butm{1}(2),'.b','MarkerSize',35)
plot(Butm{2}(1),Butm{2}(2),'.m','MarkerSize',35)
plot(Butm{3}(1),Butm{3}(2),'.g','MarkerSize',35)

% shoreline 
plot([ShoreP.B.Point1.xutm ShoreP.B.Point2.xutm],...
    [ShoreP.B.Point1.yutm ShoreP.B.Point2.yutm],'-m','LineWidth',3)
legend('Obs Direc','','Snell Direc','','','B01','B03','B05','Shoreline',...
    'location','SouthWest','fontsize',16)
title({'Swell Band: Observed vs Snells Directions',''},'fontsize',22)
ax = gca(F);
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
xlabel('UTM Easting (m)','fontsize',20)
ylabel('UTM Northing (m)','fontsize',20)
ylim([6.45*10^5,6.455*10^5])
xlim([1.7339*10^6 1.7355*10^6])


F = figure(2);clf;

quiver(Butm{2}(1),Butm{2}(2),mag*x,mag*Obs_LRy{1},'r','LineWidth',2,'MarkerSize',10)
hold on; grid on;
quiver(Butm{2}(1),Butm{2}(2),mag*x,mag*Snell_LRy{1},'k','LineWidth',2,'MarkerSize',10)
quiver(Butm{3}(1),Butm{3}(2),mag*x,mag*Obs_LRy{2},'r','LineWidth',2,'MarkerSize',10)
quiver(Butm{3}(1),Butm{3}(2),mag*x,mag*Snell_LRy{2},'k','LineWidth',2,'MarkerSize',10)

plot(Butm{1}(1),Butm{1}(2),'.b','MarkerSize',35)
plot(Butm{2}(1),Butm{2}(2),'.m','MarkerSize',35)
plot(Butm{3}(1),Butm{3}(2),'.g','MarkerSize',35)

% shoreline 
plot([ShoreP.B.Point1.xutm ShoreP.B.Point2.xutm],...
    [ShoreP.B.Point1.yutm ShoreP.B.Point2.yutm],'-m','LineWidth',3)
legend('Obs Direc','Snell Direc','','','B01','B03','B05','Shoreline',...
    'location','SouthWest','fontsize',11)
title('Sea Band: Observed vs Snells Directions','fontsize',25)
ax = gca(F);
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
xlabel('UTM Easting','fontsize',20)
ylabel('UTM Northing','fontsize',20)
ylim([6.45*10^5,6.455*10^5])


