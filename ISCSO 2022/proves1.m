clear; clc; close all;

run("DataSections.m");
x = 1:37;
x = x/37;
S = Si(:,1)/max(Si(:,1));
Smax = max(Si(:,1));
figure()
plot(x,S);
% fits: last coefficient is the term independent of x
hold on;
coef3 = polyfit(x,S,3);
y3 = polyval(coef3,x);
plot(x,y3)
plot(x,x.^3)
