%%
%%Author: Wenqiang Wang@Sustech
%%Date:2018/10/20

clear;
clc;

fid = fopen('bathy.out.linux', 'rb');
fid2 = fopen('height.dat', 'wb');
a = fread(fid, [1000, 800], 'float32');
a = fliplr(a); a = a';
fwrite(fid2, a, 'float32');
