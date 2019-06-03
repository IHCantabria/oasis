clc
clear all
close all

% path = 'Z:\Alvaro\MooringsCPP\cases\0_Clean_Case\output\';
% path = 'Z:\Alvaro\MooringsCPP\cases\0_Clean_Case_test\output\';
path = [pwd '\output\'];

nt = size(load([path 'DOF_1_Body_1.txt']),1);
[Xp,Yp] = meshgrid(-100:50:100,-100:50:100);
Zs = Xp*0.0;

B1 = load([path 'DOF_1_Body_1.txt']);
% B2 = load([path 'DOF_2_Body_1.txt']);
% B3 = load([path 'DOF_3_Body_1.txt']);
% B4 = load([path 'DOF_4_Body_1.txt']);
% B5 = load([path 'DOF_5_Body_1.txt']);
% B6 = load([path 'DOF_6_Body_1.txt']);
% 
% BB1 = load([path 'DOF_1_Body_2.txt']);
% BB2 = load([path 'DOF_2_Body_2.txt']);
% BB3 = load([path 'DOF_3_Body_2.txt']);
% BB4 = load([path 'DOF_4_Body_2.txt']);
% BB5 = load([path 'DOF_5_Body_2.txt']);
% BB6 = load([path 'DOF_6_Body_2.txt']);

X1 = load([path 'NodePosX_1.txt']);
Y1 = load([path 'NodePosY_1.txt']);
Z1 = load([path 'NodePosZ_1.txt']);
% X2 = load([path 'NodePosX_2.txt']);
% Y2 = load([path 'NodePosY_2.txt']);
% Z2 = load([path 'NodePosZ_2.txt']);
% X3 = load([path 'NodePosX_3.txt']);
% Y3 = load([path 'NodePosY_3.txt']);
% Z3 = load([path 'NodePosZ_3.txt']);
% X4 = load([path 'NodePosX_4.txt']);
% Y4 = load([path 'NodePosY_4.txt']);
% Z4 = load([path 'NodePosZ_4.txt']);

dim = [50 50 10];


figure('units','normalized','outerposition',[0 0 1 1])
hold on
grid on
box on
for it = 1:100:nt
    
    surf(Xp,Yp,Zs,'EdgeColor','none','FaceColor',[0 0.467 0.76],'FaceAlpha',0.25)

    plot3(X1(it,2:end),Y1(it,2:end),Z1(it,2:end),'ko-','linewidth',1.5,'markerfacecolor','k')
    %plot3(X2(it,2:end),Y2(it,2:end),Z2(it,2:end),'k-','linewidth',1.5)
    %plot3(X3(it,2:end),Y3(it,2:end),Z3(it,2:end),'k-','linewidth',1.5)
    %plot3(X4(it,2:end),Y4(it,2:end),Z4(it,2:end),'b-','linewidth',2.0)
    
%     plotBox(dim,[B1(it,2) B2(it,2) B3(it,2) (180/pi)*B4(it,2) (180/pi)*B5(it,2) (180/pi)*B6(it,2)],'r')
%     plotBox(dim,[BB1(it,2) BB2(it,2) BB3(it,2) (180/pi)*BB4(it,2) (180/pi)*BB5(it,2) (180/pi)*BB6(it,2)],'r')
    
    title(['Time = ' num2str(B1(it,1)) ' s'])
%     xlim([-500 500])
%     ylim([-500 500])
%     zlim([-80 20])
%     pbaspect([10 10 1])
%     xlim([-40 80])
%     ylim([-30 30])
%     zlim([-30 30])
%     pbaspect([2 1 1])
    view([45 20])
    pause(0.1)
    if it<nt
       cla
    end
end


