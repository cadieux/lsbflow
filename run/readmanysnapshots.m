function [t,dt,nx,ny,nzp,u,v,w,temp,p] = readmanysnapshots(files)
% filename = 'qavg0000.dat';
N = size(files,1);

[t,dt,nx,ny,nzp,up,vp,wp,tempp,pp] = readsnapshot(files(1).name);
u = up; v = vp; w = wp; temp = tempp; p = pp;

for i = 2:N
    [t,dt,nx,ny,nzp,up,vp,wp,tempp,pp] = readsnapshot(files(i).name);
    u = u + up;
    v = v + vp;
    w = w + wp;
    temp = temp + tempp;
    p = p + pp;
end

u = u./double(N);
v = v./double(N);
w = w./double(N);
temp = temp./double(N);
p = p./double(N);
end