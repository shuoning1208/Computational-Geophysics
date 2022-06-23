clear;clc;
tic
%% 初始化参数
nx = 154;
ny = 154;
dx = 10;
dy = 10;
nt = 2000;
dt = 0.001;
v = 3000;
fmain = 20;
e = 200;%(Gauss)
x0 = nx/2;
y0 = x0;
a = dt/dx;
X = (1:nx)*dx;
Y = (1:ny)*dy;
epoch = 50;
num_save = nt/epoch;
count = 1;

t = ((1:nt)-30)*dt;
ft = (2*(pi*fmain*t).^2-1).*exp(-(pi*fmain*t).^2);
Unew_p = zeros(nx,ny);
Unew_vx = zeros(nx,ny);
Unew_vy = zeros(nx,ny);
U_p = zeros(nx,ny);
U_vx = zeros(nx,ny);
U_vy = zeros(nx,ny);
U_p_save = zeros(nx,ny,num_save);
U_vx_save = zeros(nx,ny,num_save);
U_vy_save = zeros(nx,ny,num_save);
for it=1:nt
    U_p(x0,y0) = ft(it);
    U_vx(x0,y0) = ft(it);
    U_vy(x0,y0) = ft(it);
    for II=3:nx-2
        for JJ=3:ny-2
            dU1 = [U_p(II,JJ+1);U_vx(II,JJ+1);U_vy(II,JJ+1)]+[U_p(II,JJ-1);U_vx(II,JJ-1);U_vy(II,JJ-1)]...
                -2*[U_p(II,JJ);U_vx(II,JJ);U_vy(II,JJ)];
            dU2 = [U_p(II+1,JJ);U_vx(II+1,JJ);U_vy(II+1,JJ)]+[U_p(II-1,JJ);U_vx(II-1,JJ);U_vy(II-1,JJ)]...
                -2*[U_p(II,JJ);U_vx(II,JJ);U_vy(II,JJ)];
            coeff = -a*(-v/2-a/2*v^2);
            Um = [U_p(II,JJ);U_vx(II,JJ);U_vy(II,JJ)]+coeff*dU1+coeff*dU2;
            Unew_p(II,JJ) = Um(1,1);
            Unew_vx(II,JJ) = Um(2,1);
            Unew_vy(II,JJ) = Um(3,1); 
        end
        %边界条件
        i = 1:nx;
        j = 1:ny;
        %左边界
        Unew_p(1,j) = Unew_p(3,j);
        Unew_p(2,j) = Unew_p(3,j);
        Unew_vx(1,j) = Unew_vx(3,j);
        Unew_vx(2,j) = Unew_vx(3,j);
        Unew_vy(1,j) = Unew_vy(3,j);
        Unew_vy(2,j) = Unew_vy(3,j);
        %右边界
        Unew_p(nx,j) = Unew_p(nx-2,j);
        Unew_p(nx-1,j) = Unew_p(nx-2,j);
        Unew_vx(nx,j) = Unew_vx(nx-2,j);
        Unew_vx(nx-1,j) = Unew_vx(nx-2,j);
        Unew_vy(nx,j) = Unew_vy(nx-2,j);
        Unew_vy(nx-1,j) = Unew_vy(nx-2,j);
        %上边界
        Unew_p(i,1) = Unew_p(i,3);
        Unew_p(i,2) = Unew_p(i,3);
        Unew_vx(i,1) = Unew_vx(i,3);
        Unew_vx(i,2) = Unew_vx(i,3);
        Unew_vy(i,1) = Unew_vy(i,3);
        Unew_vy(i,2) = Unew_vy(i,3);
        %下边界
        Unew_p(i,ny) = Unew_p(i,ny-2);
        Unew_p(i,ny-1) = Unew_p(i,ny-2);
        Unew_vx(i,ny) = Unew_vx(i,ny-2);
        Unew_vx(i,ny-1) = Unew_vx(i,ny-2);
        Unew_vy(i,ny) = Unew_vy(i,ny-2);
        Unew_vy(i,ny-1) = Unew_vy(i,ny-2);
    end

    if(mod(it,100)==0)
        fprintf('time step=%d total=%d\n',it,nt);     
    end
    if(mod(it,epoch)==0)
        U_p_save(:,:,count) = Unew_p;
        U_vx_save(:,:,count) = Unew_p;
        U_vy_save(:,:,count) = Unew_p;
        count = count+1;
    end
    U_p = Unew_p;
    U_vx = Unew_vx;
    U_vy = Unew_vy;
end

pic_num = 1;
filename = 'FVM_2D.gif';
for i=1:num_save
    pcolor(X,Y,U_p_save(:,:,i))
    str = sprintf('FVM-2D\ntime=%.3f s',i*50*dt);
    xlabel('X');
    ylabel('Y');
    title(str);
    legend;
    shading interp;
    axis tight;
    axis equal;
    F = getframe(gcf);
    I = frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,filename,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
    pic_num = pic_num + 1;
end

%压力场p
figure(2)
pcolor(X,Y,U_p_save(:,:,num_save))
str = sprintf('FVM-2D\ntime=%.3f s',num_save*50*dt);
xlabel('X');
ylabel('Y');
title(str);
legend;
shading interp
axis tight;
axis equal;
xlabel('X');
ylabel('Y');
title(str)

toc
    
