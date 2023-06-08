x = 0:0.01:1;
f3 = x.^3;
f5 = x.^5;

figure(2);clf;
plot(x, f3)
hold on
plot(x, f5)