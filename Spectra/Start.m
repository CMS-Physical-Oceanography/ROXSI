%% Getting Started
%   Noah Clark
%   5/16/2023

clc;clear;

% Note: 
%   All analysis done for X01 (Asilomar; sensor 1) besides part 1
%       Asilomar - the section along the peninsula
%       sensor 1 - the sensor furthest from the shore

%% Still to do

%       - are the datetime vectors all correct and perfectly matching up? (I don't think so)
%           - I think that the wind datetime conversion converted 8 hours
%           in the wrong direction (should also be 7 instead of 8)
%           - plot Hsig to see
%       - try to make code more reusable (? ask how I could do this)
%       - break up Start.m into smaller m files
%           - upload to github
%           - comment in to describe all variables (Variables.mat)(stopped at Gfreq)
%           - remember to edit name for files where I plot; example:
%             plot_Spectro.m and load_WBdata.m
%
%
%


%% Questions I have


%       - has tide been removed from data?
%           - yes becasuse my calculated Hsigs are the same as given Hsigs?
%           - If the answer is yes or no, could any of my calculations have
%               been 
%       - how could I make the script better in general (formatting,
%       comments, etc.)?


%% Important Notes

%   - The majority of the Analysis in this script is done for X01
%   - Asilomar - the more northern section of peninsula ("X" wave buoys)
%   - China Rock - the more southern section of peninsula ("B" wave buoys)
%   - sensor 1 - the sensor furthest from the shore




%% Load in Data and Assign variables

addpath('C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)\Start Data') %folder to read in data

    %China Rock (B) Spotters:
Bspotters = {'roxsi_spotter_L2_B01_1150_reduced.mat' ...
    'roxsi_spotter_L2_B01_1158_reduced.mat' ...
    'roxsi_spotter_L2_B03_1152_reduced.mat' ...
    'roxsi_spotter_L2_B05_1153_reduced.mat'};
for i = 1:length(Bspotters)
    load(string(Bspotters{i}));
    Btime{i} = spotterL2.dtime;
    Bfreq{i} = spotterL2.frequency;
    BSee{i} = spotterL2.See;
    BGivenHsig{i} = spotterL2.Hsig;
    Bdepth{i} = -spotterL2.bottomdepth;
    BEMEM{i} = spotterL2.EMEM.meandir_f;
end

BSee{2} = cat(2,BSee{2},BSee{1});   %Buoy 1 had to be taken out of water and put back in at a later time so we must concatenate the data
BGivenHsig{2} = cat(1,BGivenHsig{2},BGivenHsig{1});
Btime{2} = cat(1,Btime{2},Btime{1});
BEMEM{2} = cat(2,BEMEM{2},BEMEM{1});
for i = 1:3     %shifting everything down (forget about {4} becuase it is now the same as {3})
    BSee{i} = BSee{i+1};
    Bfreq{i} = Bfreq{i+1};
    Btime{i} = Btime{i+1};
    BGivenHsig{i} = BGivenHsig{i+1};
    Bdepth{i} = Bdepth{i+1};
    BEMEM{i} = BEMEM{i+1};
end


    %Asilomar (X) Spotters:
Xspotters = {'roxsi_spotter_L2_X01_1151_reduced.mat' ...
    'roxsi_spotter_L2_X03_1157_reduced.mat' ...
    'roxsi_spotter_L2_X04_1155_reduced.mat'};
for i = 1:length(Xspotters)
    load(string(Xspotters{i}));
    Xtime{i} = spotterL2.dtime;
    Xfreq{i} = spotterL2.frequency;
    XSee{i} = spotterL2.See;
    XGivenHsig{i} = spotterL2.Hsig;
    Xdepth{i} = -spotterL2.bottomdepth;
    XEMEM{i} = spotterL2.EMEM.meandir_f;
end


    %General Location of China Rock and Asilomar (not loc of buoys) (taken from Google)
ChinaLAT = 36.60661;
ChinaLON = -121.9589;
AsilomarLAT = 36.6191;
AsilomarLON = -121.9376;

    %Loading Buoy Locations
load('mooringtable_largescale_2022.mat')
mooringLat = mooringtable.latitude;
mooringLon = mooringtable.longitude;

    %China Rock (B) buoy locations:
BlocIND = [25 28 31];
for i = 1:3
    Blat(i) = mooringLat(BlocIND(i));
    Blon(i) = mooringLon(BlocIND(i));
end

    %Asilomar (X) buoy locations:
XlocIND = [1 3 4];
for i = 1:3
    Xlat(i) = mooringLat(XlocIND(i));
    Xlon(i) = mooringLon(XlocIND(i));
end


    %NOAA Buoy Tide Data (Station: 9413450)(NAVD;LST;meters):
NOAA = readtable('CO-OPS_9413450_met.csv');
NOAA = NOAA(13:844,:);

    
    %Wind Data from Area:
load('roxsi_ISPAR_wind.mat')
WindDT = isparL1.dtime;
WindDT = datetime(WindDT,'TimeZone','America/Los_Angeles');
WindDT = datetime(WindDT,'TimeZone','local');
WindDir = isparL1.winddir;
WindSpd = isparL1.windspd;

    %Wind Data from Offshore:
        %this is used because the original wind data (just above) starts
        %about a week after the buoy data begins
OffWind = load('MBARI2022.mat');
%   -indecies: 3757:4537 (15 June 20:00:00 GMT  -  18 July 15:00:00 GMT)
%   -between 4537 and 4538 the data skips 10 days
OffWindTime = OffWind.NOAA.time(3757:4537);
    OffWindTime = datetime(OffWindTime,'convertfrom','datenum','TimeZone','GMT');
    OffWindTime = datetime(OffWindTime','TimeZone','America/Los_Angeles');
    OffWindTime = datetime(OffWindTime','TimeZone','local');
OffWindDir = OffWind.NOAA.wind.WDIR(3757:4537);
OffWindSpd = OffWind.NOAA.wind.WSPD(3757:4537);

    %Model See data to add to animation:
ModelSee = load('Dir_Spec_X01_SWAN.mat');
        %To save memory (in GitHub) make the Snn variable equal 0 (the
        %important variable used is Snn_f)
            ModelSee.B01.Snn = 0;
            ModelSee.X01.Snn = 0;

    %The names of all the buoys:
BBuoyName = {'B01' 'B03' 'B05'};
XBuoyName = {'X01' 'X03' 'X04'};


    %Other Variables:
df = 0.0098;
plotcolors = 'rbgmkc';  %to index when plotting in for loop



%% Spectrograms for Every Buoy

LB1x1 = [datetime('02-Jul-2022 00:00:00','TimeZone','America/Los_Angeles')...
    , datetime('02-Jul-2022 00:00:00','TimeZone','America/Los_Angeles')];
LB1x2 = [datetime('06-Jul-2022 12:00:00','TimeZone','America/Los_Angeles')...
    , datetime('06-Jul-2022 12:00:00','TimeZone','America/Los_Angeles')];
LB1y1 = [0 0.55];
LB1y2 = [0 0.55];

Gfreq = (1:129).*0.01;

        %China Rock Buoys:
        for i = 1:3
            figure(i);clf;
            delete(findall(gcf,'type','annotation'))
            set(gcf,'position',[0,570,450,450])
            pcolor(Btime{i},Gfreq,BSee{i})
            shading flat
            cb = colorbar;
            caxis([0 3])
            ylim([0 0.4])
            if i == 1
                hold on
                plot(LB1x1,LB1y1,'-r','LineWidth',2)
                plot(LB1x2,LB1y2,'-r','LineWidth',2)
                title('B01 Spectrogram')
                an = annotation('textbox',[0.4177,0.7,0.08,0.187],'string',...
                    'No Data for this Time','fontsize',8);
                set(an,'color','r')
            elseif i == 2
                title('B03 Spectrogram')
            else
                title('B05 Spectrogram')
            end
            xlabel('Date/Time')
            ylabel(cb,'Wave Energy (m^2/Hz)')
            ylabel('Frequency (Hz)')
        end
        %Asilomar Buoys:
        for i = 1:3
            figure(i + 3);clf;
            set(gcf,'position',[0,570,450,450])
            pcolor(Xtime{i},Gfreq,XSee{i})
            shading flat
            cb = colorbar;
            caxis([0 3])
            ylim([0 0.4])
            if i == 1
                title('X01 Spectrogram')
            elseif i == 2
                title('X03 Spectrogram')
            else
                title('X04 Spectrogram')
            end
            xlabel('Date/Time')
            ylabel(cb,'Wave Energy (m^2/Hz)')
            ylabel('Frequency (Hz)')
        end
        
  
        
        
%% Determining the Normal (to shoreline) Wave Direction
%Determined the shoreline points using Google Earth
    % use the formula dot(A,B) = |A||B|cos(theta)
    
ShoreP.B.Point1.Lat = 36.602756;
ShoreP.B.Point1.Lon = -121.961542;
ShoreP.B.Point2.Lat = 36.608361;
ShoreP.B.Point2.Lon = -121.958869;

ShoreP.X.Point1.Lat = 36.621764;
ShoreP.X.Point1.Lon = -121.942375;
ShoreP.X.Point2.Lat = 36.626033;
ShoreP.X.Point2.Lon = -121.940219;


proj = projcrs(26944);  %utm zone: NAD83/California zone 4
[ShoreP.B.Point1.xutm,ShoreP.B.Point1.yutm] = projfwd(proj,...
    ShoreP.B.Point1.Lat,ShoreP.B.Point1.Lon);
[ShoreP.B.Point2.xutm,ShoreP.B.Point2.yutm] = projfwd(proj,...
    ShoreP.B.Point2.Lat,ShoreP.B.Point2.Lon);
[ShoreP.X.Point1.xutm,ShoreP.X.Point1.yutm] = projfwd(proj,...
    ShoreP.X.Point1.Lat,ShoreP.X.Point1.Lon);
[ShoreP.X.Point2.xutm,ShoreP.X.Point2.yutm] = projfwd(proj,...
    ShoreP.X.Point2.Lat,ShoreP.X.Point2.Lon);

BShoreVec = [ShoreP.B.Point2.xutm - ShoreP.B.Point1.xutm,...
    ShoreP.B.Point2.yutm - ShoreP.B.Point1.yutm];
Length_BShoreVec = sqrt(BShoreVec(1)^2 + BShoreVec(2)^2);

XShoreVec = [ShoreP.X.Point2.xutm - ShoreP.X.Point1.xutm,...
    ShoreP.X.Point2.yutm - ShoreP.X.Point1.yutm];
Length_XShoreVec = sqrt(XShoreVec(1)^2 + XShoreVec(2)^2); 

NorthVec = [0 1];
Length_NorthVec = 1;

XShoreDir = acosd(dot(NorthVec,XShoreVec)/(Length_XShoreVec*Length_NorthVec));
BShoreDir = acosd(dot(NorthVec,BShoreVec)/(Length_BShoreVec*Length_NorthVec));

XNormWaveDir = 270 + XShoreDir
BNormWaveDir = 270 + BShoreDir




%% Comparing how far from normal the average directions are

B1MeanDir = mean(BEMEM{1},'all');
B1DirDiff = BNormWaveDir - B1MeanDir;

B2MeanDir = mean(BEMEM{2},'all');
B2DirDiff = BNormWaveDir - B2MeanDir;

B3MeanDir = mean(BEMEM{3},'all');
B3DirDiff = BNormWaveDir - B3MeanDir;

X1MeanDir = mean(XEMEM{1},'all');
X1DirDiff = XNormWaveDir - X1MeanDir;

X2MeanDir = nanmean(XEMEM{2},'all');
X2DirDiff = XNormWaveDir - X2MeanDir

X3MeanDir = mean(XEMEM{3},'all');
X3DirDiff = XNormWaveDir - X3MeanDir;

Location = ["China Rock";"Asilomar"];
Buoy_1_AvgDirDiff = [B1DirDiff;X1DirDiff];
Buoy_2_AvgDirDiff = [B2DirDiff;X2DirDiff];
Buoy_3_AvgDirDiff = [B3DirDiff;X3DirDiff];

NormalDirectionDifference = table(Location,Buoy_1_AvgDirDiff,...
    Buoy_2_AvgDirDiff,Buoy_3_AvgDirDiff)


%% Plotting the Locations of the Wave Buoys

    %Plot of where China Rock and where Asilomar are located on the
    %peninsula
figure(7);clf;
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
figure(8);clf;
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
figure(9);clf;
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




%% Make an time-averaged spectrum from all buoys
% I create two figures (1&2), one for China Rock and one for Asilomar

    %Determine the average depths of each buoy
for i = 1:3
    BavgD{i} = mean(Bdepth{i});
    XavgD{i} = mean(Xdepth{i});
end

    %China Rock (B):
figure(10);clf;
set(gcf,'position',[0,40,450,450])
for i = 1:3
    Bmeanspec{i} = nanmean(BSee{i}');  
    plot(Bfreq{i},Bmeanspec{i},plotcolors(i),'LineWidth',2)
    hold on
end
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Energy (m^2/Hz)')
    title('Time-Averaged Spectrum for the China Rock Buoys')
    xlim([0,0.5])
    ylim([0,1.2])
    legend(sprintf('B01 - Avg Depth (m): %g\n',round(BavgD{1},2)),...
        sprintf('B03 - Avg Depth (m): %g\n',round(BavgD{2},2)),...
        sprintf('B05 - Avg Depth (m): %g\n',round(BavgD{3},2)))

    %Asilomar (X):
figure(11);clf;
set(gcf,'position',[450,40,450,450])
for i = 1:3
    Xmeanspec{i} = nanmean(XSee{i}');
    plot(Xfreq{i},Xmeanspec{i},plotcolors(i),'LineWidth',2)
    hold on
end
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Energy (m^2/Hz)')
    title('Time-Averaged Spectrum for the Asilomar Buoys')
    xlim([0,0.5])
    ylim([0,1.2])
    legend(sprintf('X01 - Avg Depth (m): %g\n',round(XavgD{1},2)),...
        sprintf('X03 - Avg Depth (m): %g\n',round(XavgD{2},2)),...
        sprintf('X04 - Avg Depth (m): %g\n',round(XavgD{3},2)))
    
   
    
    
%% Calculate Hsig from WB spectrums and compare with given
% This is done knowing that the area under the time-averaged spectrum is
% the significant wave height

    %China Rock (B) - for each hour:
for i = 1:3
    Bt_int{i} = sum(BSee{i}).*df; %the integration
    Bt_Hsig{i} = 4*sqrt(Bt_int{i}); %m
end
    %China Rock (B) - for each buoy:
for i = 1:3
    Bbuoy_int{i} = sum(Bmeanspec{i}).*df; %the integration
    Bbuoy_Hsig{i} = 4*sqrt(Bbuoy_int{i}); %m
end

    %Asilomar (X) - for each hour:
for i = 1:3
    Xt_int{i} = sum(XSee{i}).*df; %the integration
    Xt_Hsig{i} = 4*sqrt(Xt_int{i}); %m
end
    %Asilomar (X) - for each buoy:
for i = 1:3
    Xbuoy_int{i} = sum(Xmeanspec{i}).*df; %the integration
    Xbuoy_Hsig{i} = 4*sqrt(Xbuoy_int{i}); %m
end

    %Create table with the 6 significant wave heights:
Hsig = [Bbuoy_Hsig{1};Bbuoy_Hsig{2};Bbuoy_Hsig{3};Xbuoy_Hsig{1};...
    Xbuoy_Hsig{2};Xbuoy_Hsig{3}];
Buoy = ['B01';'B03';'B05';'X01';'X03';'X04'];
Table_Hsig = table;
Table_Hsig.Buoy = Buoy;
Table_Hsig.Hsig = Hsig

    %Create figure comparing the calculated Hsigs versus the given Hsigs
        %China Rock (B):
        
figure(12);clf;
set(gcf,'position',[900,40,450,450])
x1 = datetime('15-Jun-2022','TimeZone','America/Los_Angeles');
x2 = datetime('21-Jul-2022','TimeZone','America/Los_Angeles');

for i = 1:3
    t = Btime{i};
    subplot(3,1,i)
    %xlim([x1 x2])
    
    plot(t(1:end-1),BGivenHsig{i}(1:end-1),'-b','LineWidth',1.5)
    hold on
    plot(t,Bt_Hsig{i},'-r','LineWidth',1)
    hold on
    if i == 1
        plot(t(397:398),BGivenHsig{i}(397:398),'-w','LineWidth',1.5)
        hold on
        plot(t(397:398),Bt_Hsig{i}(397:398),'-w','LineWidth',1)
        hold on
    end
    ylim([0 4]) 
    xlim([x1,x2])
    xlabel('Time (hrs)')
    ylabel('Hsig (m)')
    title({'',...
        sprintf('Given Versus Calculated Significant Wave Height for %s'...
        ,BBuoyName{i})})
    legend('Given Hsig','Calculated Hsig','location','north','fontsize',6)
    grid on
end




%% Calculate Hsig for the model spectrum
% - used trapz here but used sum and df method above
%   - I have checked that they give the same results

for i = 1:787
    M_Hsig(i) = 4*(sqrt(trapz(ModelSee.X01.f,ModelSee.X01.Snn_f(:,i),1)));
end



%% Make an animation of the spectrum from X01

    %Change the folder to save the animation figures to
cd 'C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)\Start Figures\Youtube Animation'
addpath 'C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)'

figure(20);clf
set(gcf,'position',[1350,40,450,450])
%277:444(June 27th - July 3rd)
for xx = 1:832
%   figure(20);clf;
%    delete(findall(gcf,'type','annotation'))
%     
%     if xx < 215
%        annotation('textbox',[0.32,0.8,0.4,0.07],'string','No Wind Data for this Time')
%     end
%     if xx > 214
%        subplot(3,1,1)
%        title('Wind Direction and Speed (m/s)')
%        hidden_arrow = compass(9.5,0);
%        hidden_arrow.Color = 'none';
%        hold on
%        [u,v] = pol2cart(deg2rad(WindDir(xx - 214)),WindSpd(xx - 214));
%        c = compass(u,v,'r');
%        c.LineWidth = 3;
%        ax = ancestor(c(1),'axes');
%        th = findall(ax,'Type','Text');
%        set(th,'FontSize',8)
%        %c.fontsize = 7;
%        view(90,-90)
%         annotation('textbox',[0.7,0.83,0.27,0.1],'String',...
%             sprintf('Wind Speed (m/s): %1s',...
%             sprintf('%g\n',round(WindSpd(xx - 214),2))),'fontsize',9)
%         annotation('textbox',[0.7,0.73,0.27,0.1],'String',...
%             sprintf('Wind Direction (degrees): %1s',...
%             sprintf('%g\n',round(WindDir(xx - 214),2))),'fontsize',9)
%         title({'Wind Speed and Direction',''},'fontsize',11)
%         hold off
%     end
% 
% %Plotting:
%     subplot(3,1,[2:3])
%     plot(Xfreq{1},XSee{1}(:,xx),'k','LineWidth',1.5)
% 
%     if xx > 35 && xx < 823
%         hold on
%         yy = xx - 35;
%         plot(ModelSee.X01.f,ModelSee.X01.Snn_f(:,yy),'r','LineWidth',1.5)
% 
%         annotation('textbox',[0.55,0.49,0.34,0.06],'String',...
%             sprintf('Theoretical H_1_/_3 (m): %1s',...
%             sprintf('%g\n',round(M_Hsig(yy),2))),'fontsize',9)
%         legend('Observed','Modeled','location','southeast')
%     end
% 
%     hold on
%     xlabel('Frequency (Hz)')
%     ylabel('Energy (m^2/Hz)')
%     grid on
%     axis([0 0.5 0 3.5])
%     title({'',sprintf('Frequency Spectrum for Bouy X01 @ %1s',Xtime{1}(xx))})
%     if xx < 35 && xx > 823 
%         legend('Observed','Modeled','location','southeast')
%     end    
    
    %Calculations:
        %Determining L from peak frequency/period:
            [m,indPeakFreq] = max(XSee{1}(:,xx));
            Tp(xx) = 1/Xfreq{1}(indPeakFreq);
            fp(xx) = 1/Tp(xx);
            [pWavelength(xx),pWaveNumber(xx),pCelerity(xx)] = ...
            function_wavecalculateSI(Tp(xx),XGivenHsig{1}(xx),Xdepth{1}(xx));
        %Determining L from energy weighted mean frequency/period:
            m0 = trapz(Xfreq{1},XSee{1}(:,xx),1);
            m1 = trapz(Xfreq{1},XSee{1}(:,xx).*Xfreq{1},1);
            Tm(xx) = m0/m1;
            fm(xx) = 1/Tm(xx);
            [mWavelength(xx),mWaveNumber(xx),mCelerity(xx)] = ...
            function_wavecalculateSI(Tm(xx),XGivenHsig{1}(xx),Xdepth{1}(xx));
        %Determing peak bottom orbital velocity:
            pBOV(xx) = (XGivenHsig{1}(xx)*pi)/(Tp(xx)*...
                sinh(pWaveNumber(xx)*Xdepth{1}(xx)));
        %Determing mean bottom orbital velocity:
            mBOV(xx) = (XGivenHsig{1}(xx)*pi)/(Tm(xx)*...
                sinh(mWaveNumber(xx)*Xdepth{1}(xx)));
                
%         %Adding Annotations:
%             annotation('textbox',[0.55,0.55,0.34,0.06],'String',...
%                sprintf(' Measured H_1_/_3 (m): %1s',...
%                sprintf('%g\n',round(Xt_Hsig{1}(xx),2))),'fontsize',9)
%            annotation('textbox',[0.65,0.41,0.229,0.06],'String',...
%                sprintf('L_p (m): %1s',...
%                sprintf('%g\n',round(pWavelength(xx),2))),'fontsize',9)
%             annotation('textbox',[0.65,0.35,0.229,0.06],'String',...
%                sprintf('L_m (m): %1s',...
%                sprintf('%g\n',round(mWavelength(xx),2))),'fontsize',9)
%            annotation('textbox',[0.65,0.29,0.229,0.06],'String',...
%                sprintf('u_p (m/s): %1s',...
%                sprintf('%g\n',round(pBOV(xx),2))),'fontsize',9)
%            annotation('textbox',[0.65,0.23,0.229,0.06],'String',...
%                sprintf('u_m (m/s): %1s',...
%               sprintf('%g\n',round(mBOV(xx),2))),'fontsize',9)
%     
% %To display the progress of the for-loop:
%     if xx == 0
%         disp('X01 Animation Loop is Starting')
%     elseif xx == 83
%         disp('X01 Animation Loop is 1/10 Complete')
%     elseif xx == 166
%         disp('X01 Animation Loop is 2/10 Complete')
%     elseif xx == 249
%         disp('X01 Animation Loop is 3/10 Complete')
%     elseif xx == 332
%         disp('X01 Animation Loop is 4/10 Complete')
%     elseif xx == 415
%         disp('X01 Animation Loop is 5/10 Complete')
%     elseif xx == 500
%         disp('X01 Animation Loop is 6/10 Complete')
%     elseif xx == 580
%         disp('X01 Animation Loop is 7/10 Complete')
%     elseif xx == 664
%         disp('X01 Animation Loop is 8/10 Complete')
%     elseif xx == 747
%         disp('X01 Animation Loop is 9/10 Complete')
%     elseif xx == 832
%         disp('X01 Animation Loop is Complete')
%     end
    
%%pause(0.3);clf;    %this is for running the animation in matlab
%%saveas(figure(20),sprintf('YT1_Frame%1d.jpeg',xx))
end




%% Creating multipanel time series plots for X01
% One figure showing time series plots of Hsig, energy-weighted mean and 
%  peak frequency, mean and peak bottom orbital velocity, mean and peak 
%  wavelength, and wave orbital excursion for buoy X01.

x1 = datetime('15-Jun-2022 00:00:00','TimeZone','America/Los_Angeles');
x2 = datetime('21-Jul-2022 00:00:00','TimeZone','America/Los_Angeles');

figure(13);clf;
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

    %Wave Orbital Excursion:
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




%% 3 Panel Time-Series of the spectrogram, wind speed, and tide data from NOAA for X01

    %For the entire time of recording:
figure(14);clf;
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


    %From June 27th to July 6th:
x1_2 = datetime('27-Jun-2022 00:00:00','TimeZone','America/Los_Angeles');
x2_2 = datetime('06-Jul-2022 00:00:00','TimeZone','America/Los_Angeles');

figure(15);clf;
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




%% EMEM Spectrums for Each Buoy (Directograms)

LB1x1 = [datetime('02-Jul-2022 00:00:00','TimeZone','America/Los_Angeles')...
    , datetime('02-Jul-2022 00:00:00','TimeZone','America/Los_Angeles')];
LB1x2 = [datetime('06-Jul-2022 12:00:00','TimeZone','America/Los_Angeles')...
    , datetime('06-Jul-2022 12:00:00','TimeZone','America/Los_Angeles')];
LB1y1 = [0 0.55];
LB1y2 = [0 0.55];

%China Rock Buoys:
figure(16);clf;
set(gcf,'position',[800,300,500,700])
for i = 1:3
    subplot(3,1,i)
    pcolor(Btime{i},Gfreq,BEMEM{i})
    shading flat
    cb = colorbar;
    caxis([200 360])
    ylim([0 0.55])
    if i == 1
        hold on
        plot(LB1x1,LB1y1,'-r','LineWidth',2)
        plot(LB1x2,LB1y2,'-r','LineWidth',2)
        title({'China Rock Normal Wave Direction - 292.8^o N','','','B01 Directogram'})
    elseif i == 2
        title('B03 Directogram')
    else
        title('B05 Directogram')
        xlabel('Date/Time')
    end
    ylabel(cb,'Wave Direction (degrees)')
    ylabel('Frequency (Hz)')
end

%Asilomar Buoys:
figure(17);clf;
set(gcf,'position',[950,300,500,700])
for i = 1:3
    subplot(3,1,i)
    pcolor(Xtime{i},Gfreq,XEMEM{i})
    shading flat
    cb = colorbar;
    caxis([200 360])
    ylim([0 0.55])
    if i == 1
        title({'Asilomar Normal Wave Direction - 293.9^o N','','','X01 Directogram'})
    elseif i == 2
        title('X03 Directogram')
    else
        title('X04 Directogram')
        xlabel('Date/Time')
    end
    ylabel(cb,'Wave Direction (degrees)')
    ylabel('Frequency (Hz)')
end
        
    


%% Determine the peak and mean time-weighted frequencies, bottom orbital velocities, and wavelengths for all buoys

            
% For Asilomar Buoys:
for i = 1:3
    for xx = 1:832
        %Determining L from peak frequency/period:
            [m,indPeakFreq] = max(XSee{i}(:,xx));
            XTp{i}(xx) = 1/Xfreq{i}(indPeakFreq);
            Xfp{i}(xx) = 1/XTp{i}(xx);
            [XpWavelength{i}(xx),XpWaveNumber{i}(xx),XpCelerity{i}(xx)] = ...
            function_wavecalculateSI(XTp{i}(xx),XGivenHsig{i}(xx),Xdepth{i}(xx));
        %Determining L from energy weighted mean frequency/period:
            m0 = trapz(Xfreq{i},XSee{i}(:,xx),1);
            m1 = trapz(Xfreq{i},XSee{i}(:,xx).*Xfreq{1},1);
            XTm{i}(xx) = m0/m1;
            Xfm{i}(xx) = 1/XTm{i}(xx);
            [XmWavelength{i}(xx),XmWaveNumber{i}(xx),XmCelerity{i}(xx)] = ...
            function_wavecalculateSI(XTm{i}(xx),XGivenHsig{i}(xx),Xdepth{i}(xx));
        %Determing peak bottom orbital velocity:
            XpBOV{i}(xx) = (XGivenHsig{i}(xx)*pi)/(XTp{i}(xx)*...
                sinh(XpWaveNumber{i}(xx)*Xdepth{i}(xx)));
        %Determing mean bottom orbital velocity:
            XmBOV{i}(xx) = (XGivenHsig{i}(xx)*pi)/(XTm{i}(xx)*...
                sinh(XmWaveNumber{i}(xx)*Xdepth{i}(xx)));
    end
end



% For China Rock Buoys:
for i = 2:3
    for xx = 1:832
        %Determining L from peak frequency/period:
            [m,indPeakFreq] = max(BSee{i}(:,xx));
            BTp{i}(xx) = 1/Bfreq{i}(indPeakFreq);
            Bfp{i}(xx) = 1/BTp{i}(xx);
            [BpWavelength{i}(xx),BpWaveNumber{i}(xx),BpCelerity{i}(xx)] = ...
            function_wavecalculateSI(BTp{i}(xx),BGivenHsig{i}(xx),Bdepth{i}(xx));
        %Determining L from energy weighted mean frequency/period:
            m0 = trapz(Bfreq{i},BSee{i}(:,xx),1);
            m1 = trapz(Bfreq{i},BSee{i}(:,xx).*Bfreq{1},1);
            BTm{i}(xx) = m0/m1;
            Bfm{i}(xx) = 1/BTm{i}(xx);
            [BmWavelength{i}(xx),BmWaveNumber{i}(xx),BmCelerity{i}(xx)] = ...
            function_wavecalculateSI(BTm{i}(xx),BGivenHsig{i}(xx),Bdepth{i}(xx));
        %Determing peak bottom orbital velocity:
            BpBOV{i}(xx) = (BGivenHsig{i}(xx)*pi)/(BTp{i}(xx)*...
                sinh(BpWaveNumber{i}(xx)*Bdepth{i}(xx)));
        %Determing mean bottom orbital velocity:
            BmBOV{i}(xx) = (BGivenHsig{i}(xx)*pi)/(BTm{i}(xx)*...
                sinh(BmWaveNumber{i}(xx)*Bdepth{i}(xx)));
    end
end

% - China Rock buoy B01 doesn't work for this method because there are only
%    725 values for many of the variables (I'm not sure how to get trapz to 
%    work for this)




%% Save all figures, tables, and important variables
% update this as more figures are made using Start.m

%Save figures in the "Start Figures" Folder
cd 'C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)\Start Figures'

for i = 1:17
    saveas(figure(i),sprintf('Start.Fig%1d.jpeg',i))
end

%Reassign the cd to "Start 5-16"
cd 'C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)'

%Save variables to WBvariables.mat
save('WBvariables.mat','Bdepth','BEMEM','Bfreq','BGivenHsig','Bt_Hsig',...
    'Blat','Blon','BSee','Bmeanspec','Btime','BNormWaveDir','Xdepth',...
    'XEMEM','Xfreq','XGivenHsig','Xlat','Xlon','XSee','Xtime','XavgD',...
    'M_Hsig','ModelSee','mooringLat','mooringLon','mooringtable','NOAA',...
    'OffWindDir','OffWindSpd','OffWindTime','plotcolors','WindDir',...
    'WindDT','WindSpd','NormalDirectionDifference','EXm','EXp','Xfm','Xfp',...
    'XTm','XTp','Gfreq','Hsig','XmBOV','XpBOV','XmCelerity','XpCelerity',...
    'XmWavelength','XpWavelength','ShoreP','Table_Hsig','XNormWaveDir',...
    'ChinaLAT','ChinaLON','AsilomarLAT','AsilomarLON','Xmeanspec',...
    'Xt_Hsig','Bfm','Bfp','BTm','BTp','BmBOV','BpBOV','BmCelerity',...
    'BpCelerity','BmWavelength','BpWavelength','BavgD')


close all



%% Functions Used


% function_WAVECALCULATESI - calculates an accurate estimate of the wavelength, wave
%                   number, and celerity of a function using Newton's Method.







