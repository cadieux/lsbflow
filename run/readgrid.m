function [nn,grid] = readgrid( filenames )
%READGRID Reads grid files provided to it
%   Takes in cell-array of strings of rank 3, outputs grid points

for i=1:3
    fid=fopen(filenames{i},'r');
    c = fread(fid,1,'uint32');
    nn(i) = fread(fid,1,'uint32');
    c = fread(fid,1,'uint64');
    grid{i} = fread(fid,nn(i),'float64');
    fclose(fid);
end

end

