%% NormDir.m

    % Noah Clark
    % Created: 6/16/2023
    
    % Purpose: 
    %           - to determine the incoming normal wave directions for the
    %             shoreline at China Rock and for the shoreline at Asilomar
    %           - also to create a table comparing if the wave directions
    %             get closer to normal as the waves near the shore
    %
    
    
%%

clear;clc;

    %Load in data from WBvariables.mat:
load('WBvariables.mat')

    % Select the points parallel to the beach:
        % For China Rock:
ShoreP.B.Point1.Lat = 36.602756;
ShoreP.B.Point1.Lon = -121.961542;
ShoreP.B.Point2.Lat = 36.608361;
ShoreP.B.Point2.Lon = -121.958869;
        % For Asilomar:
ShoreP.X.Point1.Lat = 36.621764;
ShoreP.X.Point1.Lon = -121.942375;
ShoreP.X.Point2.Lat = 36.626033;
ShoreP.X.Point2.Lon = -121.940219;


%%

    % Determine incoming normal wave directions:
BNormWaveDir = function_NormDir(ShoreP.B.Point1.Lat,ShoreP.B.Point1.Lon,...
    ShoreP.B.Point2.Lat,ShoreP.B.Point2.Lon)
XNormWaveDir = function_NormDir(ShoreP.X.Point1.Lat,ShoreP.X.Point1.Lon,...
    ShoreP.X.Point2.Lat,ShoreP.X.Point2.Lon)


    % Comparing how far from normal the average directions are:
B1MeanDir = mean(BEMEM{1},'all');
B1DirDiff = BNormWaveDir - B1MeanDir;
B2MeanDir = mean(BEMEM{2},'all');
B2DirDiff = BNormWaveDir - B2MeanDir;
B3MeanDir = mean(BEMEM{3},'all');
B3DirDiff = BNormWaveDir - B3MeanDir;

X1MeanDir = mean(XEMEM{1},'all');
X1DirDiff = XNormWaveDir - X1MeanDir;
X2MeanDir = nanmean(XEMEM{2},'all');
X2DirDiff = XNormWaveDir - X2MeanDir;
X3MeanDir = mean(XEMEM{3},'all');
X3DirDiff = XNormWaveDir - X3MeanDir;

    % Create Table:
Location = ["China Rock";"Asilomar"];
Buoy_1_AvgDirDiff = [B1DirDiff;X1DirDiff];
Buoy_2_AvgDirDiff = [B2DirDiff;X2DirDiff];
Buoy_3_AvgDirDiff = [B3DirDiff;X3DirDiff];

NormalDirectionDifference = table(Location,Buoy_1_AvgDirDiff,...
    Buoy_2_AvgDirDiff,Buoy_3_AvgDirDiff)


%%

    % Saving Variables:
save('WBvariables.mat','-append','NormalDirectionDifference','XNormWaveDir',...
    'BNormWaveDir','ShoreP','XNormWaveDir','BNormWaveDir')



