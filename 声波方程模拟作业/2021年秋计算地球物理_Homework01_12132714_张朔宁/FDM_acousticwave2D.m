%二维声波方程有限差分法数值模拟
clc;
clear;
tic

%分别设置空间和时间的采样点数
nx=500; 
ny=500;
nt=1000;

dx=10;    %x方向和y方向的步长相同，用dx代替
dt=0.001;   %时间步长

c=3000;    %声波在此介质中传播速度设置为3000m/s
F=20;         %震源主频设为20HZ
A=(dt*c)^2/dx^2;%计算时的常量系数
p=zeros(nx,ny,nt);%初始化数组，提高计算效率

%计算P值的循环体
for k = 2:nt-1    
    for i = 3:nx-1  
        for j = 3:ny-1
            if i==250 && j==250
                  p(i,j,k) = (1-2*(pi*F*(k*dt-1/F))^2)*exp(-(pi*F*(k*dt-1/F))^2);%在(250m，250m)处，即在二维平面的正中央选用雷克子波作为震源
            else
                p(i,j,k+1)=A*(p(i+1,j,k)+p(i-1,j,k)+p(i,j+1,k)+p(i,j-1,k)-4*p(i,j,k))-p(i,j,k-1)+2*p(i,j,k)+dt^2*(1-2*(pi*F*(k*dt-1/F))^2)*exp(-(pi*F*(k*dt-1/F))^2);
            end
        end
    end
end

%将计算的值进行成像
for k=1:10:nt
pcolor(p(:,:,k));
shading interp;
colormap('copper'); 
axis equal;%横纵坐标设置成等长刻度
axis([0,500,0,500]);
set(gca,'Ydir','reverse'); %改变坐标轴的方向
xlabel('x/m'); 
ylabel('y/m');
title(['二维声波传播图像（不设置吸收边界条件）',newline,'t=',num2str(k*dt),'s']);
image=getframe(gcf); % 捕获画面
end
toc

