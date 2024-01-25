function [NormWaveDir] = function_NormDir(SPAlat,SPAlon,SPBlat,SPBlon)
%FUNCTION_NORMWAVEDIR determines what incoming wave angle will be normal to the
%             shore in Monterey Bay, California. This function uses the formula:
%             dot(A,B) = |A||B|cos(theta)
%
%
%Input Variables:
    % SPAlat - the latitude of Point A
    % SPAlon - the longitude of Point A
    % SPBlat - the latitude of Point B
    % SPBlon - the longitude of 
%  
%Output Variables:
    % NormDir - the incoming normal direction of the waves 
%    
%Notes:
%   - PointA should be the point furthest to the west and furthest south
%   - PointB should be the point furthest to the east and furthest north
%


proj = projcrs(26944);  %utm zone: NAD83/California zone 4
[SPAxutm,SPAyutm] = projfwd(proj,SPAlat,SPAlon);
[SPBxutm,SPByutm] = projfwd(proj,SPBlat,SPBlon);


ShoreVec = [SPBxutm - SPAxutm,SPByutm - SPAyutm];
Length_ShoreVec = sqrt(ShoreVec(1)^2 + ShoreVec(2)^2);

NorthVec = [0 1]; %this is just a vector pointing directly North (angle: 0 degrees)
Length_NorthVec = 1;

ShoreDir = acosd(dot(NorthVec,ShoreVec)/(Length_ShoreVec*Length_NorthVec));

NormWaveDir = 270 + ShoreDir;

end

