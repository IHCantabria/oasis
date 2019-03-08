clc
%close all
clear all

path='Z:\Alvaro\MooringsCPP\';
color=['k-';'g-';'b-'];

ii = 2;
dt = 10;
X=load([path 'NodePosX_' num2str(ii) '.txt']);
Z=load([path 'NodePosZ_' num2str(ii) '.txt']);
[nt,ns]=size(X);

figure(1)
for i = 1:dt:nt
    hold on
    grid on
    title(['--  Time = ' num2str(X(i,1)) ' s  --'])
    plot(X(i,2:ns),Z(i,2:ns),color(ii,:),'linewidth',1.5)
    xlim([-22 2])
    ylim([-2 1])
    pause(0.01)
    cla
end

xlim([-20 1]);
ylim([-4.5 0.5]);
grid on








