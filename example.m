% Author: M73ACat
% 2023/02/06
clc
clear
close all
fs = 5120;
time = 1;
t = 0:1/fs:time-1/fs;
x1_t = 2*(1+0.5*(sin(2*pi*t))).*sin(60*pi*t);
x2_t = sin(120*pi*t);
x3_t = 0.5*cos(10*pi*t);
sig = x1_t + x2_t + x3_t;

SGCs = SGMD(sig,fs,1,0.95,0.01);

figure();
subplot(4,1,1);
plot(sig);
xlim([0,length(sig)]);
xlim([0,fs*t(end)]);
subplot(4,1,2);
plot(x1_t);
xlim([0,fs*t(end)]);
subplot(4,1,3);
plot(x2_t);
xlim([0,fs*t(end)]);
subplot(4,1,4);
plot(x3_t);
xlim([0,fs*t(end)]);


figure();
[~,b] = size(SGCs);
for i = 1:b
    subplot(b,1,i);
    plot(SGCs(:,i));
    xlim([0,fs*t(end)]);
end