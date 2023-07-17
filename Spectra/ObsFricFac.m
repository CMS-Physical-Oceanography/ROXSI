%% ObsFricFac.m

% Noah Clark 7/13/2023
% Purpose: Calculate the observed friction factor (fw_o)

% WILL EVENTUALLY BECOME A FUNCTION

% What will need to be inputted: u_br, See, h, TED


%% Example Inputs
clc;clear;
load('WBvariables.mat');
See = XSee{1};
h = Xdepth{1}; %may only need to input the time averaged depth
kw = 0.16; a1 = 5.5; a2 = -0.2; a3 = -6.3; ff = (1:129).*0.0098;
[fw,u_br,u_b] = function_FricFac(See,ff,h,kw,a1,a2,a3,'r');
TED = function_TEdis(XSee{1},XSee{2},XEMEM{1},XEMEM{2},ff,Xutm{1},...
    Xutm{2},Xdepth{1},Xdepth{2});


%% Preliminaries
rho = 1025; %kg/m^3
loop = length(ff);
df = ff(2) - ff(1);
T = 1./ff;
omegaf = 2.*pi.*ff;
See = mean(See,2); %time average
h = mean(h); %time average
%TED = mean(TED,2);

%% Calculations
a = sqrt(2.*See.*df);
a = a';

% Calculate k (by solving linear dispersion relationship):
k = zeros(1,loop);
for j = 1:loop
     Lprev = 1; Lnew = 0; thresh = 0.01; delta = 1;
     while delta > thresh
         Lprev = Lnew;
         Lnew = (9.81*T(j)^2)/(2*pi)*tanh((2*pi*h)/Lprev);
         delta = abs(Lnew - Lprev);
     end
     k(j) = (2*pi)/Lnew;
end

u_bj = a.*omegaf./(sinh(k.*h));


%fw_o = 4.*TED./(rho.*u_br.*u_bj.^2)
%the sizes of ubr and ubj do not match 

fw_o = zeros(129,832);
for i = 1:832
    for j = 1:129
        fw_o(j,i) = 4.*TED(j,i)./(rho.*u_br(i).*u_bj(j).^2);
    end
end

m_fw_o = nanmean(fw_o,2);

figure(10);clf;
plot(ff,nanmean(fw_o,2),'LineWidth',2)

