clear;clc;
tic
nx = 800;    %number of x point
nt = 1000;   % time step
c = 3000;    %velocity
density = 2500;  %define the density
shear_modulu = c^2*density;% to definite shear modulus 
dt = 0.001;
dx = 10;
e = 200;%(Gauss)
x0 = dx*nx/2;  %source position
X = (1:nx)*dx;
fmain = 10;
%% initialization
Q = zeros(2,nx);
Qnew = zeros(2,nx);
A = [0,-shear_modulu;-1/density,0];
Q_stress = zeros(nx,nt);
Q_velocity = zeros(nx,nt);
% set the source term
t = ((1:nt)-30)*dt;
ft = (2*(pi*fmain*t).^2-1).*exp(-(pi*fmain*t).^2);
%% main calculation procedure
for i=1:nt
    Qnew(1,nx/2) = ft(i);
    Qnew(2,nx/2) = ft(i);
    Q = Qnew;
    for j=2:nx-1
        dQ1 = Q(:,j+1)-Q(:,j-1);
        dQ2 = Q(:,j-1)-2*Q(:,j)+Q(:,j+1);
        Qnew(:,j) = Q(:,j) - dt/(2*dx)*A*dQ1+dt^2/(2*dx^2)*(A*A)*dQ2;
%           dQ1 = Q(:,j)-Q(:,j-1);
%           dQ2 = Q(:,j+1)-Q(:,j);
%           Qnew(:,j) = Q(:,j)-(dt/dx)*A*(dQ1);
    end
  %absorbed boundary
    Qnew(:,1) = Qnew(:,2);
    Qnew(:,nx) = Qnew(:,nx-1);

    Q_stress(:,i) = Qnew(1,:);
    Q_velocity(:,i) = Qnew(2,:);
end

pic_num = 1;
filename = sprintf('FVM_1D_fmain=%d.gif',fmain);
for i=10:10:nt
    subplot(211)
    plot(X,Q_stress(:,i)/10^6)
    str = sprintf('FVM-1D-stress\ntime step=%d',i);
    xlabel('X');
    ylabel('Stress/MPa','FontWeight','bold');
    title(str)
    subplot(212)
    plot(X,Q_velocity(:,i))
    str = sprintf('FVM-1D-velocity\ntime step=%d',i);
    title(str)
    xlabel('X');
    ylabel('Velocity/(m/s)','FontWeight','bold');
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

figure(2)
subplot(211)
plot(X,Q_stress(:,500),'linewidth',1.5)
str = sprintf('FVM-1D-stress\ntime step=%d',500);
xlabel('X');
ylabel('Stress','FontWeight','bold');
title(str)
subplot(212)
plot(X,Q_velocity(:,500),'linewidth',1.5)
str = sprintf('FVM-1D-velocity\ntime step=%d',500);
title(str)
xlabel('X');
ylabel('Velocity','FontWeight','bold');
toc
