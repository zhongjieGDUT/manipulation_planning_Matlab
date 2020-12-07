clear;
close all;
global deg
global rad
global ur10
global collision0
global collision1
global collision2
global q0
global qf
deg = 180/pi;
rad = pi/180;
%% 建立UR10模型
mdl_ur10;

ur10.qlim=[ -180 180;
            -90 90;
            -120 120;
            -90 90;
            -95 95;
            -180 180
              ]*rad;
ur10.model3d='UR/UR10';
W=[-2 2 -2.5 1.5 0 2.2];        %绘图范围
qz=[-90 0 0 0 0 0]*rad;       %初始构型
qs=[10 10 10 10 10 10]*rad;
qn=[-126 -46.8 -43.2 3.6 93.6 -36]*rad;  %逆解迭代初值
Tn=ur10.fkine(qn); Tn=se2t(Tn);
T=transl(Tn(1,4),Tn(2,4),Tn(3,4))*trotz(-pi/2)*troty(pi)*trotz(pi/2);
%工作初始造型
T0=T;

q0=ikinep(ur10,T0,qn);
qf=q0;qf(1)=-43.2*rad;
T_goal = ur10.fkine(qf);T_goal = se2t(T_goal);

%% 建立障碍模型
desk=Box(transl(0,-1.2,0.15),[1.2 0.7 0.1],'FaceColor', [255 185 15]/255, 'EdgeColor', 'k');
% 桌面的平面是0.25m
p=@(x) 2*(x-0.5).^2+0.5;
cup=Curvilinear(transl(0,-0.8,0.25),[0.3 0.3 1],p,'FaceColor',[0 191 255]/255,'EdgeColor','k');
top=Cone(transl(0,-0.8,1.25),[0.3 0.3 0.5],'FaceColor',[0 255 127]/255,'EdgeColor','k');
start_flag = Cylinder(transl(-0.7979,-0.8290,0.25),[0.2 0.2 0.001],...
    'FaceColor',[0 0.5 0.2],'EdgeColor','k');
goal_flag = Cylinder(transl(0.7179,-0.8991,0.25),[0.2 0.2 0.001],...
    'FaceColor',[0 0.99 0],'EdgeColor','k');
collision0=CollisionModel(desk);
collision1=CollisionModel(cup);
collision2=CollisionModel(top);

% %% 绘制障碍空间
% figure(1);
% desk.plot();cup.plot();start_flag.plot();goal_flag.plot();top.plot();hold on;
% ur10.plot(q0,'trail',{'b','LineWidth',2},'workspace',[-2 2 -2 2 -2 2.2],'nowrist');
% ur10.teach();
%% 三维空间
figure(2);
axis([-1.5 1.5 -1.5 1.5 0 1.5]);
hold on;
% desk.plot() ;cup.plot();top.plot();start_flag.plot();goal_flag.plot();
% plottcp(T0,'',0.0005);
% plottcp(T_goal,'',0.0005);
ur10.plot3d(qz,'trail',{'r','LineWidth',2},'nowrist');
ur10.collisions(q0,collision0);


