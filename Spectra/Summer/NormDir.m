%% NormDir.m

    % Noah Clark
    % Created: 6/16/2023
    
    % Purpose: 
    %           - to determine the incoming normal wave directions for the
    %              shoreline at China Rock and for the shoreline at Asilomar
    %           - also to create a table comparing if the wave directions
    %              get closer to normal as the waves near the shore using
    %              the average wave direction; also determine the standard
    %              deviation at each buoy for each method and add it to the
    %              values presented in the table
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
% Use my function function_NormDir

BNormWaveDir = function_NormDir(ShoreP.B.Point1.Lat,ShoreP.B.Point1.Lon,...
    ShoreP.B.Point2.Lat,ShoreP.B.Point2.Lon)
XNormWaveDir = function_NormDir(ShoreP.X.Point1.Lat,ShoreP.X.Point1.Lon,...
    ShoreP.X.Point2.Lat,ShoreP.X.Point2.Lon)


%% Using the average wave directions:

    % Comparing how far from normal the average directions are:

    % For China Rock:
BMeanDir{1} = 360 + meanangle(meanangle(BEMEM{1}(5:41,:)));
BDirDiff{1} = BNormWaveDir - BMeanDir{1};
BMeanDir{2} = 360 + meanangle(meanangle(BEMEM{2}(5:41,:)));
BDirDiff{2} = BNormWaveDir - BMeanDir{2};
BMeanDir{3} = 360 + meanangle(meanangle(BEMEM{3}(5:41,:)));
BDirDiff{3} = BNormWaveDir - BMeanDir{3};
    
for i = 1:3     % Determining standard deviation
    BavgCOLS{i} = meanangle(BEMEM{i}(5:41,:).*(pi/180));    %convert to radians
    Bavg_sumcos{i} = (sum(cos(BavgCOLS{i})))^2;
    Bavg_sumsin{i} = (sum(sin(BavgCOLS{i})))^2;
    if i == 1
        Bavg_R{i} = sqrt(Bavg_sumcos{i} + Bavg_sumsin{i})/725;
    else
        Bavg_R{i} = sqrt(Bavg_sumcos{i} + Bavg_sumsin{i})/832;
    end
    Bavg_stdev{i} = sqrt(-2*log(Bavg_R{i})).*(180/pi);  %convert back to degrees

end


    % For Asilomar:
XMeanDir{1} = 360 + meanangle(meanangle(XEMEM{1}(5:41,:)));
XDirDiff{1} = XNormWaveDir - XMeanDir{1};
XMeanDir{2} = 360 + meanangle(meanangle(XEMEM{2}(5:41,[1:172 174:end])));
XDirDiff{2} = XNormWaveDir - XMeanDir{2};
XMeanDir{3} = 360 + meanangle(meanangle(XEMEM{3}(5:41,:)));
XDirDiff{3} = XNormWaveDir - XMeanDir{3};

for i = 1:3     % Determining standard deviation
        if i == 2
            XavgCOLS{i} = meanangle(XEMEM{i}(5:41,[1:172 174:end]).*(pi/180)); %convert to radians and find time average
            Xavg_sumcos{i} = (sum(cos(XavgCOLS{i})))^2;
            Xavg_sumsin{i} = (sum(sin(XavgCOLS{i})))^2;
            Xavg_R{i} = sqrt(Xavg_sumcos{i} + Xavg_sumsin{i})/831;
            Xavg_stdev{i} = sqrt(-2*log(Xavg_R{i})).*(180/pi);
        else
            XavgCOLS{i} = meanangle(XEMEM{i}(5:41,:).*(pi/180));
            Xavg_sumcos{i} = (sum(cos(XavgCOLS{i})))^2;
            Xavg_sumsin{i} = (sum(sin(XavgCOLS{i})))^2;
            Xavg_R{i} = sqrt(Xavg_sumcos{i} + Xavg_sumsin{i})/832;
            Xavg_stdev{i} = sqrt(-2*log(Xavg_R{i})).*(180/pi);  %convert back to degrees
        end
end


    % Create Table:
Location = ["China Rock";"Asilomar"];
for i = 1:3
    AvgDirDiff{i} = [string(sprintf('%.2f %s %.2f',BDirDiff{i},...
        char(177),Bavg_stdev{i}));string(sprintf('%.2f %s %.2f',...
        XDirDiff{i},char(177),Xavg_stdev{i}))];
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
    m0 = trapz(Bfreq{i},BSee{i}(:,xx),1);
    m1 = trapz(Bfreq{i},BSee{i}(:,xx).*BEMEM{i}(:,xx));
    BEWM_dir{i}(xx) = m1/m0;
end
BEWM_DIR{i} = 360 + meanangle(BEWM_dir{i});

for i = 2:3
    for xx = 1:832
        m0 = trapz(Bfreq{i},BSee{i}(:,xx),1);
        m1 = trapz(Bfreq{i},BSee{i}(:,xx).*BEMEM{i}(:,xx));
        BEWM_dir{i}(xx) = m1/m0;
    end
    BEWM_DIR{i} = 360 + meanangle(BEWM_dir{i});
end

for i = 1:3
    BDirDiff_EWM{i} = BNormWaveDir - BEWM_DIR{i};
end

for i = 1:3     % Determining standard deviation
    BEWM_sumcos{i} = (sum(cos(BEWM_dir{i}.*(pi/180))))^2; %converting to radians
    BEWM_sumsin{i} = (sum(sin(BEWM_dir{i}.*(pi/180))))^2;
    if i == 1
        BEWM_R{i} = sqrt(BEWM_sumcos{i} + BEWM_sumsin{i})/725;
    else
        BEWM_R{i} = sqrt(BEWM_sumcos{i} + BEWM_sumsin{i})/832;
    end
    BEWM_stdev{i} = sqrt(-2*log(BEWM_R{i})).*(180/pi); %converting back to degrees
end


    % For Asilomar: 
for i = 1:3
    for xx = [1:172 174:832]
        m0 = trapz(Xfreq{i},XSee{i}(:,xx),1);
        m1 = trapz(Xfreq{i},XSee{i}(:,xx).*XEMEM{i}(:,xx));
        XEWM_dir{i}(xx) = m1/m0;
    end
    XEWM_DIR{i} = 360 + meanangle(XEWM_dir{i}([1:172 174:832]));
end

for i = 1:3
    XDirDiff_EWM{i} = XNormWaveDir - XEWM_DIR{i};
end

for i = 1:3     % Determining standard deviation
        if i == 2           
            XEWM_sumcos{i} = (sum(cos(XEWM_dir{i}.*(pi/180))))^2;   %converting to radians
            XEWM_sumsin{i} = (sum(sin(XEWM_dir{i}.*(pi/180))))^2;
            XEWM_R{i} = sqrt(XEWM_sumcos{i} + XEWM_sumsin{i})/831;
            XEWM_stdev{i} = sqrt(-2*log(XEWM_R{i})).*(180/pi);  %convert back to degrees
        else
            XEWM_sumcos{i} = (sum(cos(XEWM_dir{i}.*(pi/180))))^2;   %converting to radians
            XEWM_sumsin{i} = (sum(sin(XEWM_dir{i}.*(pi/180))))^2;
            XEWM_R{i} = sqrt(XEWM_sumcos{i} + XEWM_sumsin{i})/832;
            XEWM_stdev{i} = sqrt(-2*log(XEWM_R{i})).*(180/pi);  %convert back to degrees
        end
end


    % Create Table:
Location = ["China Rock";"Asilomar"];
for i = 1:3
    EWMDirDiff{i} = [string(sprintf('%.2f %s %.2f',BDirDiff_EWM{i},...
        char(177),BEWM_stdev{i}));string(sprintf('%.2f %s %.2f',...
        XDirDiff_EWM{i},char(177),XEWM_stdev{i}))];
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
    BPF_IND = round(Bfp{i}(xx)/df) + 1;
    BPF_dir{i}(xx) = BEMEM{i}(BPF_IND,xx);
end
BPF_DIR{1} = 360 + meanangle(BPF_dir{1});

for i = 2:3
    for xx = 1:832
        BPF_IND = round(Bfp{i}(xx)/df) + 1;
        BPF_dir{i}(xx) = BEMEM{i}(BPF_IND,xx);
    end
        BPF_DIR{i} = 360 + meanangle(BPF_dir{i});
end

for i = 1:3
    BDirDiff_PF{i} = BNormWaveDir - BPF_DIR{i};
end

for i = 1:3     % Determining standard deviation
    BPF_sumcos{i} = (sum(cos(BPF_dir{i}.*(pi/180))))^2; %converting to radians
    BPF_sumsin{i} = (sum(sin(BPF_dir{i}.*(pi/180))))^2;
    if i == 1
        BPF_R{i} = sqrt(BPF_sumcos{i} + BPF_sumsin{i})/725;
    else
        BPF_R{i} = sqrt(BPF_sumcos{i} + BPF_sumsin{i})/832;
    end
    BPF_stdev{i} = sqrt(-2*log(BPF_R{i})).*(180/pi); %converting back to degrees
end


    % For Asilomar:
for i = 1:3
    for xx = [1:172 174:832]
        XPF_IND = round(Xfp{i}(xx)/df) + 1;
        XPF_dir{i}(xx) = XEMEM{i}(XPF_IND,xx);
    end
        XPF_DIR{i} = 360 + meanangle(XPF_dir{i}([1:172 174:832]));
end

for i = 1:3
    XDirDiff_PF{i} = XNormWaveDir - XPF_DIR{i};
end

for i = 1:3     % Determining standard deviation
        if i == 2           
            XPF_sumcos{i} = (sum(cos(XPF_dir{i}.*(pi/180))))^2;   %converting to radians
            XPF_sumsin{i} = (sum(sin(XPF_dir{i}.*(pi/180))))^2;
            XPF_R{i} = sqrt(XPF_sumcos{i} + XPF_sumsin{i})/831;
            XPF_stdev{i} = sqrt(-2*log(XPF_R{i})).*(180/pi);  %convert back to degrees
        else
            XPF_sumcos{i} = (sum(cos(XPF_dir{i}.*(pi/180))))^2;   %converting to radians
            XPF_sumsin{i} = (sum(sin(XPF_dir{i}.*(pi/180))))^2;
            XPF_R{i} = sqrt(XPF_sumcos{i} + XPF_sumsin{i})/832;
            XPF_stdev{i} = sqrt(-2*log(XPF_R{i})).*(180/pi);  %convert back to degrees
        end
end


    % Create Table:
Location = ["China Rock";"Asilomar"];
for i = 1:3
    PFDirDiff{i} = [string(sprintf('%.2f %s %.2f',BDirDiff_PF{i},...
        char(177),BPF_stdev{i}));string(sprintf('%.2f %s %.2f',...
        XDirDiff_PF{i},char(177),XPF_stdev{i}))];
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



