function [Btime,Bfreq,BSee,BGivenHsig,Bdepth,BEMEM,Xtime,Xfreq,XSee,...
    XGivenHsig,Xdepth,XEMEM] = load_WBdata(anything)
%WBDATA take the spotter/wave buoy data sets and breaks them into their
%variables

%Input Variables: (no input variables are needed for this function(it just calls data))
    % anytyhing - this is not a real variable and is just a placeholder
%
%Output Variables:
    % Btime and Xtime - 
    % Bfreq and Xfreq - 
    % BSee and XSee - 
    % BGivenHsig and XGivenHsig - 
    % Bdepth and Xdepth -
    % BEMEM and XEMEM -
%
%Notes:
    % - this function only calls the data and will call the same thing
    %   every time, no matter the input
    
    
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


end