function [nn,xpts,ypts,zpts] = readgridfile( gridfile )
%READGRID Reads grid files provided to it
%   Takes in cell-array of strings of rank 3, outputs grid points

fid=fopen(gridfile,'r');
c = fread(fid,1,'uint32');
nn(1:3) = fread(fid,3,'uint32');
c = fread(fid,1,'uint64');
junk = fread(fid,5,'float64');
c = fread(fid,1,'uint64');
xpts = fread(fid,nn(1),'float64');
c = fread(fid,1,'uint64');
ypts = fread(fid,nn(2),'float64');
c = fread(fid,1,'uint64');
zpts = fread(fid,nn(3),'float64');
fclose(fid);

end