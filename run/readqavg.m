function [t,dt,xlen,ylen,zlen,...
    nx,ny,nzp,u,v,w,temp,p] = readqavg(filename)
% filename = 'qavg0000.dat';
fid=fopen(filename,'r');
%-----------------------------------
% READING PARAMETERS
%-----------------------------------
c = fread(fid,1,'uint32');
t = fread(fid,1,'float64');
dt = fread(fid,1,'float64');
% avgtperiod = fread(fid,1,'float64');
xlen = fread(fid,1,'float64');
ylen = fread(fid,1,'float64');
zlen = fread(fid,1,'float64');
c = fread(fid,2,'uint32');
nn = fread(fid,5,'uint32');
nx = nn(1); ny = nn(2); nzp = nn(3)+1;
%-----------------------------------
% READING AVG QTIES
%-----------------------------------
for l = 1:5
    c = fread(fid,1,'uint64');
    qin{l} = fread(fid,nzp*nx*ny,'float64');
    qout{l} = reshape(qin{l}, [nx,ny,nzp]);
end
fclose(fid);

u = qout{1};
v = qout{2};
w = qout{3};
temp = qout{4};
p = qout{5};


% end