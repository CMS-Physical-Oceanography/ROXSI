
load('WBvariables.mat');
See1 = XSee{1}(1:51,189);
See2 = XSee{2}(1:51,189);
Direc1 = XEMEM{1}(1:51,189);
Direc2 = XEMEM{2}(1:51,189);
ff = ff(1:51);
h1 = Xdepth{1}(189);
h2 = Xdepth{2}(189);
utm1 = Xutm{1};
utm2 = Xutm{2};


See = XSee{1}(1:51,189);
h = Xdepth{1}(189); %may only need to input the time averaged depth
kw = 0.16; a1 = 5.5; a2 = -0.2; a3 = -6.3; ff = (1:51).*0.0098;
[fwr,u_brr,u_b] = function_FricFac(See,ff,h,kw,a1,a2,a3,'r');
[fwp,u_brp] = function_FricFac(See,ff,h,kw,a1,a2,a3,'p');
[fwm,u_brm] = function_FricFac(See,ff,h,kw,a1,a2,a3,'m');
