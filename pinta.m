clc
%close all
clear all

path='Z:\Alvaro\MooringsCPP\';
color=['ko-';'go-';'bo-'];

path0 = 'C:\Users\rodriguezlua\Desktop\LABORATORIO\';
pathF1 = '0_79\';

ii = 1;
dt = 1;
ni = 1;
X=load([path 'NodePosX_' num2str(ii) '.txt']);
Z=load([path 'NodePosZ_' num2str(ii) '.txt']);
T=load([path 'CatTen_' num2str(ii) '.txt']);
[nt,ns]=size(X);
[ntT,m]=size(T);

F1 = load([path0 pathF1 'G/LabTension.dat']);

figure(3)
hold on
grid on
normT = sqrt(T(ni:ntT,5).*T(ni:ntT,5)+T(ni:ntT,6).*T(ni:ntT,6)+T(ni:ntT,7).*T(ni:ntT,7));
plot(T(ni:ntT,1),normT,'r-','linewidth',1.5)
a=0.08;
norm_T_exp = F1(:,2);
norm_T_exp = smooth(norm_T_exp,1);
plot(F1(:,1)-a,norm_T_exp,'k-','linewidth',1.5)
xlim([1 5])


%%
figure(1)
for i = 1:dt:nt
    hold on
    grid on
    title(['--  Time = ' num2str(X(i,1)) ' s  --'])
    plot(X(i,2:ns),Z(i,2:ns),color(ii,:),'linewidth',1.5)
    xlim([-10 1])
    ylim([-5 1])
    pbaspect([2 1 1])
    pause(0.1)
    if i<nt
        cla
    end
end








