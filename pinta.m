clc
%close all
clear all

figure(1)
cla

path='Z:\Alvaro\MooringsCPP\';
color=['ro-';'go-';'bo-'];
Xa=load('E:\Ordenador_Biblioteca_Alvaro\CABLE_DYNAMICS_2\Jose_Armesto\Ejemplos\3.48\prueba\Aamo_PosicionX.dat');
Za=load('E:\Ordenador_Biblioteca_Alvaro\CABLE_DYNAMICS_2\Jose_Armesto\Ejemplos\3.48\prueba\Aamo_PosicionZ.dat');
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








