clc
%close all
clear all

path='Z:\Alvaro\MooringsCPP\';
color=['k-';'g-';'b-'];

ii = 1;
dt = 10;
ni = 1;
X=load([path 'NodePosX_' num2str(ii) '.txt']);
Z=load([path 'NodePosZ_' num2str(ii) '.txt']);
T=load([path 'CatTen_' num2str(ii) '.txt']);
[nt,ns]=size(X);

figure(2)
hold on
grid on
normT = sqrt(T(ni:nt,5).*T(ni:nt,5)+T(ni:nt,6).*T(ni:nt,6)+T(ni:nt,7).*T(ni:nt,7));
plot(T(:,1),normT)


%%
figure(1)
for i = 1:dt:nt
    hold on
    grid on
    title(['--  Time = ' num2str(X(i,1)) ' s  --'])
    plot(X(i,2:ns),Z(i,2:ns),color(ii,:),'linewidth',1.5)
    xlim([-22 2])
    ylim([-5 1])
    pause(0.01)
    cla
end

xlim([-20 1]);
ylim([-4.5 0.5]);
grid on








