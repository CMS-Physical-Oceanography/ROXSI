%% FricFac2.m

%function[f_w] = function_FricFac2(See,ff,h,kw,[a1,a2,a3])

% This will eventually become that function
%       - it's not for now just for testing

%% (this should be DELETED once function completed)
clc;clear;
load('WBvariables.mat')

%% For now, assume:
% Would be entered into equation:
See = XSee{1};
h = Xdepth{1};
df = 0.0098;
kw = 0.16;  %roughness length
ff = (1:129).*df;   %vector of frequencies contained in See
tf = (1:832);   %vector of the time contained in See (not sure if I'll need)
a1 = 5.5; a2 = -0.2; a3 = -6.3; %Nielson's constants

% Other:
omegaf = 2*pi*ff;


%% Representative Orbital Velocity Method

H = 4.*sqrt(See.*df);
T(j) = 1/f(j);






