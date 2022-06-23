clear;
clc;
% 生成非结构化网格，p、e、t为返回值
% 在点矩阵p中，第一行和第二行包含网格中点的x和y坐标。
% 在边矩阵e中，第一行和第二行包含起始点和结束点的索引，第三和第四行包含起始和结束参数值，
%第五行包含边缘段编号，第六和第七行包含左侧和右侧子域编号。
% 在三角形矩阵t中，前三行包含角点的索引，按逆时针顺序给出，第四行包含子域编号。
[p, e, t] = initmesh('squareg', 'hmax', 0.1);%单元初始化
[p,e,t] = refinemesh('squareg',p,e,t);
expand = 1000; %matlab生成的网格范围默认[-1,1]，乘一个系数把范围扩大
p = expand*p;
% pdemesh(p,e,t,'nodelabels','on','elementlabels','on');%可以显示节点编号和单元编号
num_nodes = length(p);
num_elements = length(t);

%参数设置
nt = 1000; dt = 0.001;
v = 2000;
T = (1:nt)*dt;
fmain = 40;

%设置震源
s_t = (1-2*pi^2*fmain*(T-0.2).^2).*exp(-fmain*pi^2*(T-0.2).^2);

%%组装刚度矩阵

%初始化刚度矩阵
 A = zeros(num_nodes, num_nodes);
for n = 1:num_elements
%获取一个单元局部编号和全局编号之间的对应关系
    local2global = t(1:3, n);
    vertices = p(:, local2global);
     x = vertices(1, :);
     y = vertices(2, :);
     a = 0.5*((x(2)*y(3)-x(3)*y(2))-(y(3)-y(2))*x(1)+y(1)*(x(3)-x(2)));%三角形单元面积
     a = abs(a);
     a = 1/(2*a);
     A_local = zeros(3,3);
     ph = [1,0;0,1;-1,-1];
     detJ = (x(3)-x(1))*(y(2)-y(3))-(y(3)-y(1))*(x(2)-x(3));
     J_inv = a*[y(2)-y(3) , x(3)-x(2);y(3)-y(1) , x(1)-x(3)];
            for i = 1: 3
                for j = 1: 3
                    A_local(i, j) = -0.5*detJ.*(ph(i,:)*J_inv*(J_inv'*ph(j,:)'));
                end
            end
            A(local2global, local2global) = A(local2global, local2global) + A_local;
            if(mod(n,100)==0)
                fprintf('number of  elements:%d/%d\n',n,num_elements)
            end
 end
    

%组装质量矩阵和右端向量
 M = zeros(num_nodes, num_nodes); % 总质量矩阵
 Me = zeros(3, 3);
 F = zeros(num_nodes, 1);% 初始化右端向量
 Fe = zeros(3, 1); % 初始单元右端向量

 % 单元质量矩阵
        N{1} = @(xi, eta) 1 - xi - eta; N_xi{1} = -1; N_eta{1} = -1;
        N{2} = @(xi, eta) xi; N_xi{2} = 1; N_eta{2} = 0;
        N{3} = @(xi, eta) eta; N_xi{3} = 0; N_eta{3} = 1;
        ymax = @(xi) 1 - xi;

 for i = 1: num_elements
   local2global = t(1: 3, i);%获取当前单元包含的节点编号
   vertices = p(:, local2global);%获取当前单元所有节点的x,y坐标
   xx = vertices(1, :); yy = vertices(2, :);
   J = [xx(1) * N_xi{1} + xx(2) * N_xi{2} + xx(3) * N_xi{3}, xx(1) * N_eta{1} + xx(2) * N_eta{2} + xx(3) * N_eta{3};
                yy(1) * N_xi{1} + yy(2) * N_xi{2} + yy(3) * N_xi{3}, yy(1) * N_eta{1} + yy(2) * N_eta{2} + yy(3) * N_eta{3}];
            detJ = abs(det(J));  
       Me(1,1) = 1/12*detJ;Me(1,2) = 1/24*detJ;Me(1,3) = 1/24*detJ;
            Me(2,1) = 1/24*detJ;Me(2,2) = 1/12*detJ;Me(2,3) = 1/24*detJ;
            Me(3,1) = 1/24*detJ;Me(3,2) = 1/24*detJ;Me(3,3) = 1/12*detJ;
            Fe(:,1) = 1/6*detJ;
            M(local2global, local2global) = M(local2global, local2global) + Me;
            F(local2global, 1) = F(local2global, 1) + Fe;
            if(mod(i,100)==0)
                fprintf('number of  elements:%d/%d\n',i,num_elements)
            end
 end


boundarynodes = unique([e(1, :) e(2, :)]);
U = zeros(num_nodes,nt);
[m,source_x] = find(abs(p(1,:))>=0&abs(p(1,:))<=50&abs(p(2,:))>=0&abs(p(2,:))<=50);

%利用递推关系式求波场值
for i = 2:nt-1
    U(source_x(1),i) = s_t(i);
    right_vector = v^2*dt^2*A*U(:,i);
    U(:,i+1) = 2*U(:,i)-U(:,i-1)-M\right_vector;
    U(boundarynodes,i+1) = 0;
    if(mod(i,10)==0)
        fprintf('time step=%d total=%ds\n',i,nt);
    end
end

filename = sprintf('FEM-2D dt=%.3f fm=%.1f v=%d explicit.gif',dt,fmain,v);
pic_num = 1;
for i=5:5:nt
    trisurf(t(1: 3, :)', p(1, :)', p(2, :)', U(:,i))
    str = sprintf('time step=%d\ndt=%.3f fm=%.1f v=%d',i,dt,fmain,v);
    title(str)
    colorbar
    shading interp
    view([90,90])
    F = getframe(gcf);
    I = frame2im(F,256);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,filename,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
    pic_num = pic_num + 1;
end