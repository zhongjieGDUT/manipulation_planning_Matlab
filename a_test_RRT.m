%% 用 RRT 解决避障抓取问题 IB-RRT*存在问题，所以从头开始
% Data: 2019/10/23
% author:Jie Zhong
% 调试成功！
clear;
close all;
%% 初始化
% 环境建模
setur10; %数据建立
% close all;
%%%%% parameters end here %%%%%
source = q0;        %初始位型
goal   = qf;        %目标位型
N = 6000;            %采样总次数

stepsize =    36*rad;alfa=1;steer=1;
disTh    =    36*rad;
if ~feasibleJoint(source),error('source lies on an obstacle');end
if ~feasibleJoint(goal),error('source lies on an obstacle');end
RRTree = [source -1];  %  RRT rooted at the source, representation node and parent index
pathFound = false;
h = waitbar(0,'采样中');
tic;
%% 主循环__原始版本

for i=1:N
    str = ['N=',num2str(i)];
    waitbar(i/N,h,str);
    if rand < 0.5
        x_rand = sample();
    else
        x_rand = goal; % sample taken as goal to bias tree generation to goal
    end
    
    [~, I]=min(distanceCost(RRTree(:,1:6),x_rand) ,[],1); 
    % A 是距离sample最近的节点距离，I是其索引，find closest as per the function in the metric node to the sample
    closestNode = RRTree(I,1:6);
    %%%%%%%%%%%%%% Steer 的法一
    if steer
        delta = x_rand - closestNode;
        d = delta./sqrt(sum(delta.^2));
        newPoint = closestNode+stepsize*d;
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Steer 的法二
    else
        T1=ur10.fkine(closestNode);
        T2=ur10.fkine(x_rand);
        e=tr2delta(T1,T2);
        dq=pinv(ur10.jacobe(closestNode))*e;
        newPoint = closestNode+alfa*dq';
        %%%%%%%%%%%%%%
    end
    
    t = ur10.fkine(newPoint).t';
    if t(3)<=0,continue;end
    
  
    [~, I2]=min(distanceCost(RRTree(:,1:6),newPoint),[],1);
    
    if distanceCost(newPoint,RRTree(I2,1:6))<disTh, continue; end 
    
    if ~checkPath([RRTree(I2,1:6);newPoint]) % if extension of closest node in tree to the new point is feasible
        continue;
    end
    if distanceCost(newPoint,goal)<disTh && checkPath([newPoint;goal])
       pathFound=true;
       RRTree=[RRTree;newPoint,I2]; % add node 
       break; 
    end % goal reached
    
    RRTree=[RRTree;newPoint I2]; % add node 
    t=ur10.fkine(newPoint).t';
    disp(num2str(length(RRTree)));
    plot3(t(1),t(2),t(3),'b.','MarkerSize',1.5); %画出落点
    
end
%% 主循环__版本一
% while ~pathFound
%     %采样
%     i = int32(rand*vol);
%     x_rand = X_free(i,:);
%     
%     [~, I]=min(distanceCost(RRTree(:,1:6),x_rand(1:6)) ,[],1); 
%     closestNode = RRTree(I,1:6);
%     delta = x_rand(1:6) - closestNode ;
%     d = ((delta-min(delta))/(max(delta)-min(delta))-0.5)*2;
%     newPoint = closestNode+stepsize*d;
%     
%     [A, I2]=min(distanceCost(RRTree(:,1:6),newPoint),[],1);
%     if distanceCost(newPoint,goal)<disTh, pathFound=true;break; end % goal reached
%     if ~checkPath([RRTree(I2,1:6);newPoint]) % if extension of closest node in tree to the new point is feasible
%         continue;
%     end
%     
%     
%     
%     
%     
%     
% end
%% 记录规划的
path = [goal];
prev=size(RRTree,1);
while prev>0
    path=[RRTree(prev,1:6);path];
    prev=RRTree(prev,7);
end
%% 可视化路径
figure(2);
for i=1:size(path,1)
    name=['UR1#',num2str(i)];
    ro = creat_ur(name);
    ro.plot3d(path(i,:),'noname','nowrist','workspace',W);
    if ~feasibleJoint(path(i,:))
        disp('fail');break;
    end
end
T_end=size(path,1);
t=0:T_end-1;
xx=0:0.5:T_end-1;
n=length(xx);
tr=zeros(n,6);
figure('Name','RRT joint trajectory ');hold on;grid on;
for i=1:6
    y  = path(:,i)';
    cs = spline(t,[0 y 0]);
    plot(xx,ppval(cs,xx));
    tr(:,i)=ppval(cs,xx)';
end
hold off;
figure('Name','RRT');hold on;
for i=1:size(tr,1)
    name=['UR#+',num2str(i)];
    ro = creat_ur(name);
    ro.plot3d(tr(i,:),'tile2color',[0.5 1 0.5],'noname','nowrist','workspace',W);
end
desk.plot() ;cup.plot();top.plot();start_flag.plot();goal_flag.plot();hold off;
%% 需要用的子函数
%  这里定义了几个子函数
%% 子函数1: 判断关节角是否可行
function flag=feasibleJoint(q) %输入关节角行向量
    global ur10
    global collision0
    global collision1
    global collision2
    flag = true;
    %判断机器人与环境是否有碰撞
    if ur10.collisions(q,collision0) ...
            ||  ur10.collisions(q,collision1) || ...
            ur10.collisions(q,collision2)
        flag=false;
    end
end
%% 子函数2：采样函数，返回可行范围内的可行关节角度
function x_rand=sample()
    global ur10
    qlim = ur10.qlim;  %获取关节角范围
    m = ur10.n;        %采样次数
    for i=1:m
        x_rand(i) = qlim(i,1)+rand*(qlim(i,2)-qlim(i,1));
    end
end
%% 子函数3：获得采样点的邻域节点集合
function X_near = NearestVertices(x_rand,RRTree)
% NearestVertices 
% 返回树RRTree中距离采样点距离r以内的点
% 记录其坐标，在树中的坐标，父节点。
    
    [num,~] = size(RRTree);
    gama = 2*pi;
    r = gama*(log(num)/num)^(1/3);
    dis = distanceCost(RRTree(:,1:6),x_rand);
    index = find(dis<=r);
    [n,~] = size(index);
    if n==0
        X_near = [];
    else
        X_near = zeros(n,8);
        for i=1:n
            X_near(i,:) = [RRTree(index(i),1:7),index(i)];
        end
    end 
end
%% 子函数4：计算每个节点到新的采样点的距离，返回距离向量
function h=distanceCost(a,b) %返回列向量
    h = sqrt(sum((a-b).^2,2));
end
%% 子函数5：若邻域集合为空集，则取一个最近的节点
function X_near = NearestVertex(x_rand,RRTree)
    dis = distanceCost(RRTree(:,1:6),x_rand);
    index=find(min(dis));
    X_near=[RRTree(index,1:7),index];
end
%% 子函数6：GetSortedList()
function Ls = GetSortedList(x_new,X_near,RRTree)
    [n,~]=size(X_near);
    Ls=cell(1,3,n);
    step=2;
    for i=1:n
        traj = Steer(X_near(i,1:6),x_new,step);
        cost = Cost(RRTree,X_near(i,:))+distanceCost(x_new,X_near(i,1:6));
        Ls{1,1,i} = X_near(i,:);
        Ls{1,2,i} = cost;
        Ls{1,3,i} = traj;
    end
    Ls = sortList(Ls);
end
%% 子函数7: 在两个关节角之间等距插值
function traj = Steer(x_near,x_new,step)
    traj=zeros(step,6);
    for i=1:6
        traj(:,i)=linspace(x_near(i),x_new(i),step)';
    end
end
%% 子函数8 ：计算每个节点的代价
function  pathLength = Cost(RRTree,x_near)
    %COST 计算从根节点到x_near的路径长度
    path=[x_near(1:6)];
    prev=x_near(7); 
    pathLength=0;
    if prev<2
        pathLength=pathLength+distanceCost(x_near,RRTree(1,1:6));
    else
         while prev>0
            path=[RRTree(prev,1:6);path];
            prev=RRTree(prev,7);
         end
         for i=1:size(path,1)-1
            pathLength=pathLength+distanceCost(path(i,1:6),path(i+1,1:6));
         end
    end
end
%% 子函数9：按代价给邻域集合排序
function sorted_Ls = sortList(Ls)
    [~,~,n]=size(Ls);
    sorted_Ls = cell(1,3,n);
    cost=zeros(n,1);
    for i=1:n
        cost(i)=Ls{1,2,i};
    end
    [~,index]=sort(cost);
    for i=1:n
        for j=1:3
            sorted_Ls{1,j,i}=Ls{1,j,index(i)}; 
        end
    end
    
end
%% 子函数10：选择最好的父节点
function [x_min,cost] = ChooseBestParentTree(Ls)
    x_min = [];
    cost = inf;
    [~,~,n] = size(Ls);
    for i=1:n
        traj=Ls{1,3,i};%获取路径
        if checkPath(traj)
            x_min = Ls{1,1,i};
            cost = Ls{1,2,i};
            break;
        end
    end
end
%% 子函数11：看路径上是否无碰
function feasible = checkPath(traj)
    feasible = true;
    step=15; 
    tau = zeros(step,6);
%     for i=1:6
%         tau(:,i) = linspace(traj(1,i),traj(2,i),step)';
%     end
    tau = jtraj(traj(1,:),traj(2,:),step);
    for i=1:step
        if ~feasibleJoint(tau(i,:))
            feasible=false;
            break;
        end
    end
end
%% 子函数12：RewireVertices
function RRTree_new=RewireVertices(x_min,x_new,Ls,RRTree)
    [~,~,n]=size(Ls);
    for i=1:n
        x_near=Ls{1,1,i};
        traj=Ls{1,3,i};
        dis=Cost(RRTree,x_min);
        if dis+distanceCost(x_min(1:6),x_new)+distanceCost(x_new,x_near(1:6))<Cost(RRTree,x_near)
            if checkPath([x_near(1:6);x_new])
                RRTree(x_near(8),7)=size(RRTree,1);
            end
        end
    end
    RRTree_new=RRTree;
end