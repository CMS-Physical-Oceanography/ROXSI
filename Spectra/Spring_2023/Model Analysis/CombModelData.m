%% CombModelData.m
%
% Noah Clark            5/14/2024
%
% Purpose: To load in the provided data from the 4 models (SWAN050, 
%          SWAN050sm, COUP050, and SWAN075), add them to their own
%          structures, and then save them in ModelResults.mat
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


%% Preliminaries

clc;clear;
addpath('Model Results')


%% SWAN050

load('SWAN2022_Rogkn050_sub.mat')

SWAN050.B01 = B01;
SWAN050.B03 = B03;
SWAN050.B05 = B05;
SWAN050.B10 = B10;
SWAN050.B13 = B13;
SWAN050.X01 = X01;
SWAN050.X03 = X03;
SWAN050.X04 = X04;
SWAN050.X05 = X05;
clear B01 B03 B05 B07 B10 B13 X01 X03 X04 X05 X08


%% SWAN050sm

load('SWAN2022_Madkn050_sub_sm.mat')

SWAN050sm.B01 = B01;
SWAN050sm.B03 = B03;
SWAN050sm.B05 = B05;
SWAN050sm.B10 = B10;
SWAN050sm.B13 = B13;
SWAN050sm.X01 = X01;
SWAN050sm.X03 = X03;
SWAN050sm.X04 = X04;
SWAN050sm.X05 = X05;
clear B01 B03 B05 B07 B10 B13 X01 X03 X04 X05 X08


%% COUP050

load('COUP2022_Rogkn050_sub.mat')

COUP050.B01 = B01;
COUP050.B03 = B03;
COUP050.B05 = B05;
COUP050.B10 = B10;
COUP050.B13 = B13;
COUP050.X01 = X01;
COUP050.X03 = X03;
COUP050.X04 = X04;
COUP050.X05 = X05;
clear B01 B03 B05 B07 B10 B13 X01 X03 X04 X05 X08


%% SWAN075

load('SWAN2022_Rogkn075_sub.mat')

SWAN075.B01 = B01;
SWAN075.B03 = B03;
SWAN075.B05 = B05;
SWAN075.B10 = B10;
SWAN075.B13 = B13;
SWAN075.X01 = X01;
SWAN075.X03 = X03;
SWAN075.X04 = X04;
SWAN075.X05 = X05;
clear B01 B03 B05 B07 B10 B13 X01 X03 X04 X05 X08


%% Save


save('ModelResults.mat','SWAN050','SWAN050sm','COUP050','SWAN075')








