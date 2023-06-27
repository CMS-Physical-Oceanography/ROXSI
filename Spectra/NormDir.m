%% NormDir.m

    % Noah Clark
    % Created: 6/16/2023
    
    % Purpose: 
    %           - to determine the incoming normal wave directions for the
    %              shoreline at China Rock and for the shoreline at Asilomar
    %           - also to create a table comparing if the wave directions
    %              get closer to normal as the waves near the shore using
    %              the average wave direction
    %                   - make the same table but using the mean energy 
    %                      weighted wave direction
    %                   - again, make the same table but instead using the
    %                      wave direction of the peak frequency
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


df = 0.0098;


%% Determine incoming normal wave directions:

BNormWaveDir = function_NormDir(ShoreP.B.Point1.Lat,ShoreP.B.Point1.Lon,...
    ShoreP.B.Point2.Lat,ShoreP.B.Point2.Lon)
XNormWaveDir = function_NormDir(ShoreP.X.Point1.Lat,ShoreP.X.Point1.Lon,...
    ShoreP.X.Point2.Lat,ShoreP.X.Point2.Lon)


%% Using the average wave directions:

    % Comparing how far from normal the average directions are:

BMeanDir{1} = 360 + meanangle(meanangle(BEMEM{1}(5:41,:)));
BDirDiff{1} = BNormWaveDir - BMeanDir{1};

BMeanDir{2} = 360 + meanangle(meanangle(BEMEM{2}(5:41,:)));
BDirDiff{2} = BNormWaveDir - BMeanDir{2};

BMeanDir{3} = 360 + meanangle(meanangle(BEMEM{3}(5:41,:)));
BDirDiff{3} = BNormWaveDir - BMeanDir{3};

XMeanDir{1} = 360 + meanangle(meanangle(XEMEM{1}(5:41,:)));
XDirDiff{1} = XNormWaveDir - XMeanDir{1};

XMeanDir{2} = 360 + meanangle(meanangle(XEMEM{2}(5:41,[1:172 174:end])));
XDirDiff{2} = XNormWaveDir - XMeanDir{2};

XMeanDir{3} = 360 + meanangle(meanangle(XEMEM{3}(5:41,:)));
XDirDiff{3} = XNormWaveDir - XMeanDir{3};

    % Create Table:
Location = ["China Rock";"Asilomar"];
for i = 1:3
    AvgDirDiff{i} = [BDirDiff{i};XDirDiff{i}];
end

AvgNormDirDiff = table(Location,AvgDirDiff{1},...
    AvgDirDiff{2},AvgDirDiff{3});
AvgNormDirDiff.Properties.VariableNames{2} = 'Buoy #1';
AvgNormDirDiff.Properties.VariableNames{3} = 'Buoy #2';
AvgNormDirDiff.Properties.VariableNames{4} = 'Buoy #3'


%% Using the energy weighted average wave directions:

    % For China Rock:
i = 1;
for xx = 1:725
    BEWM_IND = round(Bfm{i}(xx)/df);
    BEWM_dir{i}(xx) = BEMEM{i}(BEWM_IND,xx);
end
BEWM_DIR{1} = 360 + meanangle(BEWM_dir{1});

for i = 2:3
    for xx = 1:832
        BEWM_IND = round(Bfm{i}(xx)/df);
        BEWM_dir{i}(xx) = BEMEM{i}(BEWM_IND,xx);
    end
        BEWM_DIR{i} = 360 + meanangle(BEWM_dir{i});
end

for i = 1:3
    BDirDiff_EWM{i} = BNormWaveDir - BEWM_DIR{i};
end

    % For Asilomar:
for i = 1:3
    for xx = [1:172 174:832]
        XEWM_IND = round(Xfm{i}(xx)/df);
        XEWM_dir{i}(xx) = XEMEM{i}(XEWM_IND,xx);
    end
        XEWM_DIR{i} = 360 + meanangle(XEWM_dir{i}([1:172 174:832]));
end

for i = 1:3
    XDirDiff_EWM{i} = XNormWaveDir - XEWM_DIR{i};
end

    % Create Table:
Location = ["China Rock";"Asilomar"];
for i = 1:3
    EWMDirDiff{i} = [BDirDiff_EWM{i};XDirDiff_EWM{i}];
end

EWMNormDirDiff = table(Location,EWMDirDiff{1},...
    EWMDirDiff{2},EWMDirDiff{3});
EWMNormDirDiff.Properties.VariableNames{2} = 'Buoy #1';
EWMNormDirDiff.Properties.VariableNames{3} = 'Buoy #2';
EWMNormDirDiff.Properties.VariableNames{4} = 'Buoy #3'


%% Using the wave direction of the peak frequency:

    % For China Rock:
i = 1;
for xx = 1:725
    BPF_IND = round(Bfp{i}(xx)/df);
    BPF_dir{i}(xx) = BEMEM{i}(BPF_IND,xx);
end
BPF_DIR{1} = 360 + meanangle(BPF_dir{1});

for i = 2:3
    for xx = 1:832
        BPF_IND = round(Bfp{i}(xx)/df);
        BPF_dir{i}(xx) = BEMEM{i}(BPF_IND,xx);
    end
        BPF_DIR{i} = 360 + meanangle(BPF_dir{i});
end

for i = 1:3
    BDirDiff_PF{i} = BNormWaveDir - BPF_DIR{i};
end

    % For Asilomar:
for i = 1:3
    for xx = [1:172 174:832]
        XPF_IND = round(Xfp{i}(xx)/df);
        XPF_dir{i}(xx) = XEMEM{i}(XPF_IND,xx);
    end
        XPF_DIR{i} = 360 + meanangle(XPF_dir{i}([1:172 174:832]));
end

for i = 1:3
    XDirDiff_PF{i} = XNormWaveDir - XPF_DIR{i};
end


    % Create Table:
Location = ["China Rock";"Asilomar"];
for i = 1:3
    PFDirDiff{i} = [BDirDiff_PF{i};XDirDiff_PF{i}];
end

PFNormDirDiff = table(Location,PFDirDiff{1},...
    PFDirDiff{2},PFDirDiff{3});
PFNormDirDiff.Properties.VariableNames{2} = 'Buoy #1';
PFNormDirDiff.Properties.VariableNames{3} = 'Buoy #2';
PFNormDirDiff.Properties.VariableNames{4} = 'Buoy #3'


%% Saving Variables

save('WBvariables.mat','-append','AvgNormalDirectionDifference','XNormWaveDir',...
    'BNormWaveDir','ShoreP','XNormWaveDir','BNormWaveDir','BMeanDir',...
    'XMeanDir')



