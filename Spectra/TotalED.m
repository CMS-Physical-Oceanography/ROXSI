%% TotalED.m


    % Noah Clark
    % Created: 6/21/2023
    
    % Purpose: 
    %           - calculate the values used to determine the total energy
    %              dissipation and then calculate/estimate the total energy
    %              dissipation due to wave breaking and bottom friction
    %           - determine the flux through each buoy at each frequency 
    %              and compare to the flux (at same frequency) that went
    %              through the previous buoy (previous --> more offshore)
    %           - determine the distance and angle between each of the
    %              neighboring buoys
    %           - determine delta_r for each buoy 
    
    
%%

clear;clc;

    %Load in data from WBvariables.mat:
load('WBvariables.mat')



%% Determine Flux
% EX: BFluxDiff{1} represents the difference in flux from buoy B01 to buoy
%     B03 (positive number means a decrease in flux)

    % Constants:
df = 0.0098;
rho = 1025; %kg/m^3
g = 9.81; %m/s^2


    % For China Rock Buoys:
for i = 1:3
    h = BavgD{i};
    for j = 1:129
        f = j*df;
        T = 1/f;
        omega = 2*pi*f;

        H = 4*sqrt(var(BSee{i}(j,:)));
        a = H/2;
        BE{i}(j) = 0.5*rho*g*a^2;
        [L,k,WDP,WS,RD,C] = function_wavecalculateSI(T,H,h);

        BCg{i}(j) = 0.5*(1 + (2*k*h)/(sinh(2*k*h)))*(omega/k);

        BFlux{i}(j) = BE{i}(j)*BCg{i}(j);
    end
end
BFluxDiff{1} = BFlux{1} - BFlux{2};
BFluxDiff{2} = BFlux{2} - BFlux{3};


    % For Asilomar Buoys:
for i = 1:3
    h = XavgD{i};  
    for j = 1:129
        f = j*df;
        T = 1/f;
        omega = 2*pi*f;
        if i == 2
            H = 4*sqrt(var(XSee{i}(j,[1:172 174:end])));
        else 
            H = 4*sqrt(var(XSee{i}(j,:)));
        end
        a = H/2;
        XE{i}(j) = 0.5*rho*g*a^2;
        [L,k,WDP,WS,RD,C] = function_wavecalculateSI(T,H,h);

        XCg{i}(j) = 0.5*(1 + (2*k*h)/(sinh(2*k*h)))*(omega/k);

        XFlux{i}(j) = XE{i}(j)*XCg{i}(j);
    end
    
XFluxDiff{i} = diff(XFlux{i});
end
XFluxDiff{1} = XFlux{1} - XFlux{2};
XFluxDiff{2} = XFlux{2} - XFlux{3};



%% Determine delta_r (delta_r = L*cos(theta))

% 1) take the given x and y coordinates from the original spotter buoy data
%     set
% 2) determine the distance between with pythagorean theorem
% 3) determine wave propagation direction between two buoys (take average
%     of the 2 average directions)
% 4) determine theta by comparing the wave propagation direction to the
% angle formed between the line connecting the two sites
% 5) plug all into formula to find delta_r


    % Add the path for the data:
addpath('C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)\Start Data')


            % For China Rock Buoys:
    Bspotters = {'roxsi_spotter_L2_B01_1150_reduced.mat' ...
        'roxsi_spotter_L2_B03_1152_reduced.mat' ...
        'roxsi_spotter_L2_B05_1153_reduced.mat'};
for i = 1:3
    load(string(Bspotters{i}));
    Bxx{i} = spotterL2.X;
    Byy{i} = spotterL2.Y;
end

    %Determining the distance between buoys:
%  **EX: BDist(1) = distance between buoys 1 and 2**
BDist(1) = sqrt((Bxx{2} - Bxx{1})^2 + (Byy{2} - Byy{1})^2);
BDist(2) = sqrt((Bxx{3} - Bxx{2})^2 + (Byy{3} - Byy{2})^2);

    % Converting to utm coordinates:
proj = projcrs(26944);
[B01utmx,B01utmy] = projfwd(proj,Blat(1),Blon(1));
[B03utmx,B03utmy] = projfwd(proj,Blat(2),Blon(2));
[B05utmx,B05utmy] = projfwd(proj,Blat(3),Blon(3));


    % Determining the direction of one buoy to the next (relative to
    % North):
BV1x = B03utmx - B01utmx;
BV1y = B03utmy - B01utmy;
BV1 = [BV1x BV1y];
Length_BV1 = sqrt(BV1x^2 + BV1y^2);

BV2x = B05utmx - B03utmx;
BV2y = B05utmy - B03utmy;
BV2 = [BV2x BV2y];
Length_BV2 = sqrt(BV2x^2 + BV2y^2);

NorthVec = [0 1];
Length_NorthVec = 1;

dir_BV1 = acosd(dot(NorthVec,BV1)/(Length_BV1*Length_NorthVec));
dir_BV2 = acosd(dot(NorthVec,BV2)/(Length_BV2*Length_NorthVec));

for j = 1:129
    BMDir{1}(j) = meanangle(BEMEM{1}(j,:));
    BMDir{2}(j) = meanangle(BEMEM{2}(j,:));
    BWaveDir{1}(j) = mean([BMDir{1}(j) BMDir{2}(j)]) - 180;
    Btheta{1}(j) = BWaveDir{1}(j) - dir_BV1 + 360; 
    
    BMDir{3}(j) = meanangle(BEMEM{3}(j,:));
    BWaveDir{2}(j) = mean([BMDir{1}(j) BMDir{2}(j)]) - 180;     %the wave direction is assumed to be the average of the wave directions at each of the two buoys
    Btheta{2}(j) = BWaveDir{2}(j) - dir_BV2 + 360;
end
Bdelta_r{1} = BDist(1).*cosd(Btheta{1});
Bdelta_r{2} = BDist(2).*cosd(Btheta{2});



            % For Asilomar Buoys:
    Xspotters = {'roxsi_spotter_L2_X01_1151_reduced.mat' ...
        'roxsi_spotter_L2_X03_1157_reduced.mat' ...
        'roxsi_spotter_L2_X04_1155_reduced.mat'};
for i = 1:3
    load(string(Xspotters{i}));
    Xxx{i} = spotterL2.X;
    Xyy{i} = spotterL2.Y;
end

    % Determining the distance between buoys:
XDist(1) = sqrt((Xxx{2} - Xxx{1})^2 + (Xyy{2} - Xyy{1})^2);
XDist(2) = sqrt((Xxx{3} - Xxx{2})^2 + (Xyy{3} - Xyy{2})^2);

    % Converting to utm coordinates:
proj = projcrs(26944);
[X01utmx,X01utmy] = projfwd(proj,Xlat(1),Xlon(1));
[X03utmx,X03utmy] = projfwd(proj,Xlat(2),Xlon(2));
[X04utmx,X04utmy] = projfwd(proj,Xlat(3),Xlon(3));

    % Determining the direction of one buoy to the next (relative to
    % North):
XV1x = X03utmx - X01utmx;
XV1y = X03utmy - X01utmy;
XV1 = [XV1x XV1y];
Length_XV1 = sqrt(XV1x^2 + XV1y^2);

XV2x = X04utmx - X03utmx;
XV2y = X04utmy - X03utmy;
XV2 = [XV2x XV2y];
Length_XV2 = sqrt(XV2x^2 + XV2y^2);

NorthVec = [0 1];
Length_NorthVec = 1;

dir_XV1 = acosd(dot(NorthVec,XV1)/(Length_XV1*Length_NorthVec));
dir_XV2 = acosd(dot(NorthVec,XV2)/(Length_XV2*Length_NorthVec));

for j = 1:129
    XMDir{1}(j) = meanangle(XEMEM{1}(j,:));
    XMDir{2}(j) = meanangle(XEMEM{2}(j,[1:172 174:end]));
    XWaveDir{1}(j) = mean([XMDir{1}(j) XMDir{2}(j)]) - 180;
    Xtheta{1}(j) = XWaveDir{1}(j) - dir_XV1 + 360; 
    
    XMDir{3}(j) = meanangle(XEMEM{3}(j,:));
    XWaveDir{2}(j) = mean([XMDir{1}(j) XMDir{2}(j)]) - 180;     %the wave direction is assumed to be the average of the wave directions at each of the two buoys
    Xtheta{2}(j) = XWaveDir{2}(j) - dir_XV2 + 360;
end
Xdelta_r{1} = XDist(1).*cosd(Xtheta{1});
Xdelta_r{2} = XDist(2).*cosd(Xtheta{2});



%% Estimate the spatially averaged rate of wave energy dissipation
% - This energy dissipation is due to both wave breaking and bottom friction
% - the term accounting for any decrease in wave energy flux resulting from
%    two-dimensional wave refraction effects (delta_R_j) isn't used to
%     calculate this at the moment


    % For China Rock Buoys:
B_TED{1} = BFluxDiff{1}./Bdelta_r{1};
B_TED{2} = BFluxDiff{2}./Bdelta_r{2};


    % For Asilomar Buoys:
X_TED{1} = XFluxDiff{1}./Xdelta_r{1};
X_TED{2} = XFluxDiff{2}./Xdelta_r{2};




%% Saving Variables

B01utm = [B01utmx B01utmy]; Butm{1} = B01utm;
B03utm = [B03utmx B03utmy]; Butm{2} = B03utm;
B05utm = [B05utmx B05utmy]; Butm{3} = B05utm;
X01utm = [X01utmx X01utmy]; Xutm{1} = X01utm;
X03utm = [X03utmx X03utmy]; Xutm{2} = X03utm;
X04utm = [X04utmx X04utmy]; Xutm{3} = X04utm;


save('WBvariables.mat','-append','B_TED','X_TED','Butm','Xutm')


    



