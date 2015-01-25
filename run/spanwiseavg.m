function [u,v,w,temp,p] = spanwiseavg(nx,ny,nzp,uin,vin,win,tempin,pin)
% computes spanwise average
u(1:nx,1:nzp)=sum(uin,2)./double(ny); 
v(1:nx,1:nzp)=sum(vin,2)./double(ny);
w(1:nx,1:nzp)=sum(win,2)./double(ny); 
temp(1:nx,1:nzp)=sum(tempin,2)./double(ny);
p(1:nx,1:nzp)=sum(pin,2)./double(ny);
end