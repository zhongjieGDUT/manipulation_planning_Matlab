%% 用bi_RRT* 解决避障抓取问题
% Data: 2019/10/30
% author:Jie Zhong
% 优化的目标（代价函数）是整个过程关节变化量最小

clear;
close all;
%% 初始化
% 环境建模
setur10; %数据建立
%%%%% parameters end here %%%%%
source = q0;        %初始位型
goal   = qf;        %目标位型
stepsize = 36*rad;  %步长
N = 1000;            %采样总次数
pl = zeros(1,N);    %记录每次迭代所找到的路径长度，形成路径进化曲线
if ~feasibleJoint(source),error('source lies on an obstacle');end
if ~feasibleJoint(goal),error('source lies on an obstacle');end
RRTree1 = [source -1];  % First RRT rooted at the source, representation node and parent index
RRTree2 = [goal -1];    % Second RRT rooted at the goal, representation node and parent index
pathFound=false;
cost_best = inf;
h = waitbar(0,'采样中');
tic;
%% 主循环
for i=1:N
    pl(i) = cost_best;
    flag1=0;flag2=0;
    str = ['N=',num2str(i)];
    waitbar(i/N,h,str);
    %% 采样
    x_rand = sample();
    %%  对 RRTree1 操作
    [~, I]=min(distanceCost(RRTree1(:,1:6),x_rand) ,[],1);
    x_nearest= RRTree1(I,1:6);
    delta = x_rand - x_nearest;
    d = delta./sqrt(sum(delta.^2));
    % 获取采样点x_new
    x_new = x_nearest + stepsize.*d;
    if ~feasibleJoint(x_new),continue;end
    
    X_near = NearestVertices(x_new,RRTree1);
    if isempty(X_near)
        X_near = NearestVertex(x_new,RRTree1);
    end
    Ls = GetSortedList(x_new,X_near,RRTree1);
    x_min = ChooseBestParent(Ls);

    if ~isempty(x_min)
        RRTree1=[RRTree1;[x_new,x_min(8)]];
        RRTree1=RewireVertices(x_min,x_new,Ls,RRTree1);
        flag1=1;
    end
    [~, I1]=min(distanceCost(RRTree2(:,1:6),x_new) ,[],1);
    x_connect=RRTree2(I1,1:6);
    
    [flag2,RRTree2]=ConnectGraphs(RRTree2,x_new,x_connect,stepsize);
    if flag1==1 && flag2==1
        pathlen = Cost(RRTree1,RRTree1(end,:))+Cost(RRTree2,RRTree2(end,:));
        if pathlen<cost_best
            cost_best = pathlen;
            disp(['find in step=',num2str(i)]);
            tree1=RRTree1;tree2=RRTree2;pl(i)=cost_best;
        end
    end
    a=RRTree1;
    b=RRTree2;
    RRTree1=b;
    RRTree2=a;
end
%% 绘制路径计算路径长度
path=[tree1(end,1:6)];
prev=tree1(end,7);
while prev>0
    path=[tree1(prev,1:6);path];
    prev=tree1(prev,7);
end
prev=tree2(end,7);
path1=[tree2(prev,1:6)];
prev=tree2(prev,7);
while prev>0
    path1=[tree2(prev,1:6);path1];
    prev=tree2(prev,7);
end  
path=[path;flipud(path1)];
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
xx=0:0.1:T_end-1;
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
    x_rand = zeros(1,6); 
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
    gama = pi*2;
    r = gama*(log(num)/num)^(1/2);        %这个参数值得商榷
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
    
    h =sqrt(sum((a-b).^2,2));
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
        pathLength=pathLength+distanceCost(x_near(1:6),RRTree(1,1:6));
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
function x_min = ChooseBestParent(Ls)
    x_min = [];
    [~,~,n] = size(Ls);
    for i=1:n
        traj=Ls{1,3,i};%获取路径
        if checkPath(traj)
            x_min = Ls{1,1,i};
            break;
        end
    end
end
%% 子函数11：看路径上是否无碰
function feasible = checkPath(traj)
    feasible = true;
    step=10; 
    tau = zeros(step,6);
    for i=1:6
        tau(:,i) = linspace(traj(1,i),traj(2,i),step)';
    end
%     tau = jtraj(traj(1,:),traj(2,:),step);
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
%% 子函数13：ConnectGraphs
function [flag,RRTree]=ConnectGraphs(RRTree,x1,x2,stepsize)
    delta = x1-x2;
    d = delta./sqrt(sum(delta.^2));
    x_new = x2+stepsize*d;
    X_near=NearestVertices(x_new,RRTree);
    Ls = GetSortedList(x1,X_near,RRTree);
    x_min = ChooseBestParent(Ls);
%     traj=[];
    flag=0;
    if ~isempty(x_min)
        RRTree = [RRTree;[x1,x_min(8)]];
        flag=1;
    end
end