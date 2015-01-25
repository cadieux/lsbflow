function [t,dt,nx,ny,nzp,u,v,w,temp,p] = readsnapshot(filename)
% clear all;
% filename = 'q00001.dat';
fid=fopen(filename,'r');
%-----------------------------------
% READING PARAMETERS
%-----------------------------------
prebuf(1:7) = fread(fid,7,'float32');
% c = fread(fid,1,'uint32');
% t = fread(fid,1,'float32');
% dt = fread(fid,1,'float32');
% avgtperiod = fread(fid,1,'float32');
% c = fread(fid,2,'uint32');
% nn = fread(fid,3,'uint32');
% nx = nn(1); ny = nn(2); nzp = nn(3)+1;
isum = int32(prebuf(1));
t = prebuf(2);
dt = prebuf(3);
nx = int32(prebuf(4));
ny = int32(prebuf(5));
nzp = int32(prebuf(6));
nypl = int32(prebuf(7));
% c = fread(fid,1,'uint32');
%-----------------------------------
% READING AVG QTIES
%-----------------------------------
% for l = 1:5
res = nx*nzp*ny;
qin = fread(fid,5*res,'float32');
qout = reshape(qin,[nx,nzp,ny,5]);
clear qin;

% u = reshape(qout(:,:,:,1),[nx,nypl,nzp]);
% ii = 1;
for j = 1: ny
    for k = 1: nzp
        for i = 1: nx
            u(i,j,k) = qout(i,k,j,1);
            v(i,j,k) = qout(i,k,j,2);
            w(i,j,k) = qout(i,k,j,3);
            temp(i,j,k) = qout(i,k,j,4);
            p(i,j,k) = qout(i,k,j,5);
        end
    end
end
clear qout;

%     qout{l} = reshape(qin{l},[nx,nzp,ny]);
% end
fclose(fid);

% for j=1:ny
%     u(:,j,:) = qout{1}(:,:,j);
%     v(:,j,:) = qout{2}(:,:,j);
%     w(:,j,:) = qout{3}(:,:,j);
%     temp(:,j,:) = qout{4}(:,:,j);
%     p(:,j,:) = qout{5}(:,:,j);
% end

% end