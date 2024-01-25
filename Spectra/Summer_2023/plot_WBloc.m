%% plot_WBloc.m

    % Noah Clark
    % Created: 6/16/2023
    
    % Purpose: plot the general locations of China Rock and Asilomar and
    %          the locations of the wave buoys 
 
    
    
%%

clear;clc;

    %Load in data from WBvariables.mat:
load('WBvariables.mat')

    %General Location of China Rock and Asilomar (not loc of buoys)(taken
    % from Google):
ChinaLAT = 36.60661;
ChinaLON = -121.9589;
AsilomarLAT = 36.6191;
AsilomarLON = -121.9376;


%%

    %Plot of where China Rock and where Asilomar are located on the
    %peninsula
figure(1);clf;
set(gcf,'position',[450,570,450,450])
% nlat = [36.6 36.63];
%  nlon = [-121.96 -121.93];
%  geolimits(nlat,nlon)
geoscatter(ChinaLAT,ChinaLON,'sb','linewidth',10)
hold on
geoscatter(AsilomarLAT,AsilomarLON,'sr','linewidth',10);
geobasemap satellite
legend('China Rock','Asilomar','location','northwest')
title('General Location of China Rock and Asilomar')
 

    %Plot of location of the China Rock Buoys
figure(2);clf;
set(gcf,'position',[900,570,450,450])
for i = 1:3
    geoscatter(Blat(i),Blon(i),'o',plotcolors(i),'linewidth',5)
    hold on
end
        %Adding the shoreline points
        geoplot([ShoreP.B.Point1.Lat,ShoreP.B.Point2.Lat],...
            [ShoreP.B.Point1.Lon,ShoreP.B.Point2.Lon],'m','linewidth',2)
        geoscatter([ShoreP.B.Point1.Lat,ShoreP.B.Point2.Lat],...
            [ShoreP.B.Point1.Lon,ShoreP.B.Point2.Lon],'k','linewidth',3)
geobasemap satellite
title('Location of the China Rock Wave Buoys')
legend({'B01','B03','B05',['Parallel' newline 'Shoreline']},...
    'orientation','horizontal','fontsize',7)


    %Plot of location of the Asilomar Buoys
figure(3);clf;
set(gcf,'position',[1350,570,450,450])
for i = 1:3
    geoscatter(Xlat(i),Xlon(i),'o',plotcolors(i),'linewidth',5)
    hold on
end
        %Adding the shoreline points
        geoplot([ShoreP.X.Point1.Lat,ShoreP.X.Point2.Lat],...
            [ShoreP.X.Point1.Lon,ShoreP.X.Point2.Lon],'m','linewidth',2)
        geoscatter([ShoreP.X.Point1.Lat,ShoreP.X.Point2.Lat],...
            [ShoreP.X.Point1.Lon,ShoreP.X.Point2.Lon],'k','linewidth',3)
geobasemap satellite
title('Location of the Asilomar Wave Buoys')
legend({'X01','X03','X04',['Parallel' newline 'Shoreline']},...
    'orientation','horizontal','fontsize',7)


%%

    % Saving Variables:
save('WBvariables.mat','-append','ChinaLAT','ChinaLON','AsilomarLAT','AsilomarLON')




