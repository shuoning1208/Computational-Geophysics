clear
clc
tic

dt = 0.001; % 单位（s）
ds = 10; % 空间差分：ds = dx = dz 
nx = 800;%800个空间采样点
nz = 800;
nt = 2000;%2000个时间采样点 故总长为2s
fmain = 20;%设置震源的主频为20HZ

slice = 10;%每间隔10步保存结果
slice_num = nt/slice;

t = 1:nt;%计算震源（生成一个数组）
t0 = 50;%设置偏移量
s_t = (1-2*(pi*fmain*dt*(t-t0)).^2).*exp(-(pi*fmain*dt*(t-t0)).^2);
 
%计算过程中数组变量初始化
P_current = zeros(nx,nz);
P_next = zeros(nx,nz);
P_past = zeros(nx,nz);
Px_fft = complex(zeros(nx,nz));
Pz_fft = complex(zeros(nx,nz));
Px_ifft = complex(zeros(nx,nz));
Pz_ifft = complex(zeros(nx,nz));
P = zeros(nx,nz);
P_slice =zeros(nx,nz,slice_num);

fid = fopen('bathy.out.linux', 'rb');
fid2 = fopen('height.dat', 'wb');
a = fread(fid, [800, 800], 'float32');
a = fliplr(a); a = a';%高程文件
fwrite(fid2, a, 'float32');
V = sqrt(abs(a)*9.8);

slice_count = 1;

%计算波数
kmax = pi/ds;
delta_kx = kmax/(nx/2);
delta_kz = kmax/(nz/2);
kx = zeros(nx,1);
kz = zeros(nz,1);
i = 1:nx/2;
j = 1:nz/2;
kx(i,1) = i*delta_kx;
kx(nx/2+i,1) = -kmax + i*delta_kx;
kz(j,1) = j*delta_kz;
kz(nz/2+j,1) = -kmax + j*delta_kz;

for T = 1:nt
    P_current(400,400) = P_current(400,400) + dt^2*s_t(1,T);%将震源的位置设置在坐标（400，400）

    %傅里叶变换
    Px_fft = fft(P_current');%这里主要要进行转置的运算
    Pz_fft = fft(P_current);

    i = 1:nz;
    Px_fft(:,i) = (-kx.^2).*Px_fft(:,i);%(ik)^2*p(k)
    i = 1:nx;
    Pz_fft(:,i) = (-kz.^2).*Pz_fft(:,i);

    %反傅里叶变换
    Px_ifft = ifft(Px_fft)';
    Pz_ifft = ifft(Pz_fft);
    P = real(Px_ifft+Pz_ifft);%取实部

    %递推公式计算
    P_next = 2*P_current-P_past+(V.^2*dt^2).*P;
    P_past = P_current;
    P_current = P_next;

    if(mod(T,slice)==0)
        P_slice(:,:,slice_count) = P_next;%将采样切片数据存储
        slice_count = slice_count + 1;
    end
end
%绘制图像
filename = '海啸波.gif';
figure(1)
for i = 1:slice_num
    imagesc(P_slice(:,:,i));
    shading interp;
    axis equal;
    colormap('hot');
    set(gca,'yDir','reverse');
    str_title = ['海啸波  t=',num2str(dt*i*10),'s'];
    title(str_title);
    colorbar;
    drawnow; %刷新屏幕
    F = getframe(gcf);%捕获图窗作为影片帧
    I = frame2im(F); %返回图像数据
    [I, map] = rgb2ind(I, 256); %将rgb转换成索引图像
    if i == 1
        imwrite(I,map, filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(I,map, filename,'gif','WriteMode','append','DelayTime',0.1);
    end
    
end
toc