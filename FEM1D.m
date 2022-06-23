clear;
clc;
%网格大小、数量设置
dx = 10;
dt = 0.001; 
nx = 500; 
nt = 1000;    
v = 2000; 
fmain = 40;
%变量初始化设置
x = (1:nx)*dx;
t = (1:nt)*dt;
u = zeros(nx,nt);
%加载震源，雷克子波
st = (1-2*pi^2*fmain*(t-0.5).^2).*exp(-fmain*pi^2*(t-0.5).^2);
u(floor(nx/2),:) = st;
%设置边界值
u(1,:) = 0;
u(nx,:) = 0;

%初始化刚度矩阵和质量矩阵
stiffness_matrix = zeros(nx-2,nx-2);
mass_matrix = zeros(nx-2,nx-2);

F = zeros(nx-2,1);
%计算刚度矩阵和质量矩阵
for i=1:nx-2
    if i>1
        stiffness_matrix(i-1,i) = -1/dx;
        mass_matrix(i-1,i) = dx/6;
        stiffness_matrix(i,i-1) = -1/dx;
        mass_matrix(i,i-1) = dx/6;   
    end
    stiffness_matrix(i,i) = 2/dx;
    mass_matrix(i,i) = 2/3*dx;
end

%计算时间项
for j = 2:nt-1
    u(floor(nx/2),j) = st(j);
    u(1,j) = 0;
    u(nx,j) = 0;
%         load_vector =  mass_matrix * (2*u(2:end-1,j)-u(2:end-1,j-1));
%     u(2:end-1,j+1) =  ( dt^2*v^2 * stiffness_matrix + mass_matrix ) \ ( load_vector );
    right = v^2*dt^2*stiffness_matrix*u(2:end-1,j);
    u(2:end-1,j+1) = 2*u(2:end-1,j)-u(2:end-1,j-1)-mass_matrix\right;
end

figure(1)
n = size(u,2);
num_slice = size(10:10:n,1);
fmat=moviein(num_slice);
filename =sprintf('FEM_1D_fmain=%d.gif',fmain);
pic_num = 1;
for i=10:10:n
    plot(x,u(:,i));
    axis([-inf inf -1 1]);
    grid on;
    str = sprintf('time step=%d,dt=%.3f,dx=%d,fmain=%.1f',i,dt,dx,fmain);
    title(str);
    fmat(:,i)=getframe;
    f=getframe(gcf);
    I=frame2im(f);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,filename,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
    pic_num = pic_num + 1;
end

figure(2)
plot(x,u(:,nt));
str = sprintf('time step=%d,dx=%d dt=%.3f v=%d fm=%d',nt,dx,dt,v,fmain);
title(str);
saveas(gcf, 'FEM-1D', 'png')