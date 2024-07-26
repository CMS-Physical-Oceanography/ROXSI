%% calc_EFractionM.m
%
% Noah Clark            5/3/2024
%
%
% Purpose: To calculate the percentage of energy in the sea and swell bands
%          at each of the observation points in the model. These
%          percentages are averaged to determine what percentage of the
%          energy is in the sea and swell bands on average.
%
% Note:
% - The frequency indecies for Sea: 0.088 - 0.252 Hz       (ff(13)-ff(29))
%       - for observed Sea: 0.088 - 0.245 Hz
% - The frequency indecies for Swell: 0.040 - 0.088 Hz     (ff(1)-ff(13))
%       - for observed Swell: 0.0392 - 0.088 Hz
%
% - the total fraction adds above 100% (101.7%) because they share the
%   frequency 0.088 Hz
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


%% Preliminaries

clc;clear;

addpath('Model Results')
load('ModelResults.mat')

ff = Model.B01.f';


%% Should I also make the DD plot for B03 like I did for obs?
% ???
% ???
% ???


%% Calculate the average percentages

% B01
B01FracSea = sum(Model.B01.WS(13:29,:),'all')/sum(Model.B01.WS,'all'); % 0.835
B01FracSwell = sum(Model.B01.WS(1:13,:),'all')/sum(Model.B01.WS,'all'); % 0.177
% B03
B03FracSea = sum(Model.B03.WS(13:29,:),'all')/sum(Model.B03.WS,'all'); % 0.800
B03FracSwell = sum(Model.B03.WS(1:13,:),'all')/sum(Model.B03.WS,'all'); % 0.219
% B05
B05FracSea = sum(Model.B05.WS(13:29,:),'all')/sum(Model.B05.WS,'all'); % 0.744
B05FracSwell = sum(Model.B05.WS(1:13,:),'all')/sum(Model.B05.WS,'all'); % 0.274

% B10
B10FracSea = sum(Model.B10.WS(13:29,:),'all')/sum(Model.B10.WS,'all'); % 0.738
B10FracSwell = sum(Model.B10.WS(1:13,:),'all')/sum(Model.B10.WS,'all'); % 0.281
% B13
B13FracSea = sum(Model.B13.WS(13:29,:),'all')/sum(Model.B13.WS,'all'); % 0.722
B13FracSwell = sum(Model.B13.WS(1:13,:),'all')/sum(Model.B13.WS,'all'); % 0.294

% Average for all 5
AVGFracSea = mean([B01FracSea B03FracSea B05FracSea B10FracSea B13FracSea]); % 0.768
AVGFracSwell = mean([B01FracSwell B03FracSwell B05FracSwell B10FracSwell B13FracSwell]); % 0.249


%% Create Table

Buoys = ['B01';'B03';'B05';'B10';'B13';'AVG'];
SeaFraction = [B01FracSea;B03FracSea;B05FracSea;B10FracSea;B13FracSea;AVGFracSea];
SeaFraction = round(SeaFraction,3);
SwellFraction = [B01FracSwell;B03FracSwell;B05FracSwell;B10FracSwell;B13FracSwell;AVGFracSwell];
SwellFraction = round(SwellFraction,3);

MSeaSwellFraction = table(Buoys,SeaFraction,SwellFraction)









