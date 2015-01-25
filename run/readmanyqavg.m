function [t,dt,xlen,ylen,zlen,...
    nx,ny,nzp,u,v,w,temp,p] = readmanyqavg(files)
% filename = 'qavg0000.dat';
N = size(files,1);
for i = 1:N
    fid=fopen(files(i).name,'r');
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
    
    if (i==1)
        u = qout{1};
        v = qout{2};
        w = qout{3};
        temp = qout{4};
        p = qout{5};
    else
        u = u + qout{1};
        v = v + qout{2};
        w = w + qout{3};
        temp = temp + qout{4};
        p = p + qout{5};
    end
end

u = u./N;
v = v./N;
w = w./N;
temp = temp./N;
p = p./N;

end