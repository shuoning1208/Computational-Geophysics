%设计一个函数对给定不同的参数都可以生成相应的二维声波响应
tic
p = acoustic2D(20,2000,800,800,2000,400,400,10,10,0.001);
toc
%函数主体
function p1 = acoustic2D( f, nt, nx, nz, c, sx, sz, dx, dz, dt )
    %对变量进行初始化，提高计算的效率
    source = zeros(nt,1);
    pnew = zeros(nx,nz);
    pold = zeros(nx,nz);
    p = zeros(nx, nz);
    dx2 = zeros(nx, nz);
    dz2 = zeros(nx, nz);
    
    for i = 1:nt
            source(i)=(1-(2*pi*f*(i*dt-1/f))^2)*exp(-(pi*f*(i*dt-1/f))^2);  %随时间的推移，源的幅值会发生变化   
        for k= 2:nz-1
            for j = 2:nx-1
                dx2(j,k) = (p(j+1,k)-2*p(j,k)+p(j-1,k))/(dx^2);%space derivative
                dz2(j,k) = (p(j,k+1)-2*p(j,k)+p(j,k-1))/(dz^2);
                pnew(j,k) = c^2*(dt^2)*(dx2(j,k)+dz2(j,k))+2*p(j,k)-pold(j,k);%迭代计算
            end
        end 
        pnew(sx,sz) = pnew(sx,sz)+(dt^2)*source(i);% add source term : sx means the x position of source,the same as sz
        pold = p;%将现在的P值赋给下一循环的pold
        p = pnew;%将pnew的值赋给下一循环的p
     
        imagesc(p);%对p这一二维矩阵中相应元素的值按照颜色深浅显示出来(实时)
        axis equal;
        title([  'FD for 2D acousic wave propagation', newline,'t = ',num2str(i*dt),'s']);
        xlabel('X distance/m');
        ylabel('Z distance/m'); 
        drawnow  
        p1 = p;
    end
end
   
  
  
    
            
            