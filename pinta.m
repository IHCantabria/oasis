clc
%close all
clear all

figure(1)
cla

path='/home/biblioteca/GIT_repositorios/NuevoFEM_cpp/';
color=['ro-';'go-';'bo-'];
Xa=load('/home/biblioteca/Desktop/CABLE_DYNAMICS/Jose_Armesto/Ejemplos/3.48/prueba/Aamo_PosicionX.dat');
Za=load('/home/biblioteca/Desktop/CABLE_DYNAMICS/Jose_Armesto/Ejemplos/3.48/prueba/Aamo_PosicionZ.dat');
nn=size(Xa,2);
for ii=1:3
    X=load([path 'NodePosX_' num2str(ii) '.txt']);
    Z=load([path 'NodePosZ_' num2str(ii) '.txt']);

    n=size(X,2);

    plot(X(2:n),Z(2:n),color(ii,:))
    hold on
    clear X
    clear Z
end
%plot(Xa(1,2:nn),Za(1,2:nn),'ko--')
xlim([-20 1]);
ylim([-4.5 0.5]);
grid on








