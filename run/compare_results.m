% plot and compare results between postproc and otfavg
clc;
clear all;
% close all;

% qavgfile = 'qavg00012.dat';
qavgfiles = dir('q00029*');
% otffile = 'otfproc0003.dat';
gridfile = 'grid.dat';
% gridfiles = {'xgrid.dat';'ygrid.dat';'zgrid.dat'};

% read in grid
% [nn,grid] = readgrid(gridfiles);
% xpts = grid{1};
% ypts = grid{2};
% zpts = grid{3};
[dim,xpts,ypts,zpts] = readgridfile(gridfile);

% read in otf avg qties
% [t,nx,ny,nzp,cfo,utauo,cfbl,cpo,...
%     dstaro,thetao,dstarb,thetab,upo,zpo] = readotfproc(otffile);

% read in time avg qties
% [t,dt,xlen,ylen,zlen,...
%     nx,ny,nzp,uin,vin,win,tempin,pin] = readqavg(qavgfile);

% [t,dt,nx,ny,nzp,uin,vin,win,tempin,pin] = readsnapshot(qavgfile);
[t,dt,nx,ny,nzp,uin,vin,win,tempin,pin] = readmanysnapshots(qavgfiles);

% ubar(1:nzp) = uin(1,1,:);
% 
% for k=1:nzp
%     xi(k) = zpts(1)*0.5*( 1 + cos(pi*double(k-1)/double(nzp-1)) );
%     z3
% end
% 
% ubar2 = interp1(zpts,ubar,z2,'spline');

% [t,dt,xlen,ylen,zlen,...
%     nx,ny,nzp,uin,vin,win,tempin,pin] = readmanyqavg(qavgfiles);

% for i=1:nx
%     up(i,:,:) = uin(i,:,:) - uin(nx-3,:,:);
% end
% 
% pin = pin/dt - 0.5*(up.^2 + vin.^2 + win.^2);


% load Spalart data for comparisons
load('SpalartDNS.mat');

% perform spanwise avg
[u,v,w,temp,p] = spanwiseavg(nx,ny,nzp,uin,vin,win,tempin,pin);

% missing p factor
% p = 11/6.*p;

% some parameters
nz = nzp-1;
xnu = 0.001;
u0 = 1.;
den = 0.5*u0^2;

% if (size(xpts)~=nx)
%     for i=1:nx
%         xpts(i) = 10^5/6*xnu + 10.5*10^5/3.*xnu*(i-1)/nx;
%     end
% end

% calc Cf
% dz1 = zpts(nzp) - zpts(nzp-1);
% dz2 = zpts(nzp-1) - zpts(nzp-2);
% dudz_w(1:nx) = (4*u(:,nzp-1) - u(:,nzp-2) -3*u(:,nzp))./(dz1+dz2);
% utau = xnu*dudz_w;

dz = zpts(nz)-zpts(nzp);
utau = xnu*(u(:,nz)-u(:,nzp))/dz;

cf = utau/den;
% cf1 = utau(:,1)./(0.5*u(:,1).^2);

dzplus = zpts(nz)*sqrt(0.5*abs(cf))/xnu;
% calc Cp
% first normalize pressure
% for k=1:nzp
%     p0(k) = abs(min(p(:,k)));
%     pnorm(:,k) = p(:,k) + p0(k);
% end
% p0 = abs(min(min(p)));
% pnorm = p0 + p;
% pnorm = pnorm./max(max(pnorm));
% cp = (pnorm(:,nzp) - pnorm(:,1))/den;
%  -min(p(:,nz-17))
% p_ind = find(zpts/zpts(1)>=0.1,1,'last');
p_ind = nzp;
% cp = (p(:,p_ind) - p(1,nx/10))/den;
cp = (p(:,p_ind) - min(p(1:nx/4,p_ind)))/den;
% cp = (p(:,p_ind) - p(10,1))/den;
% cp = cp - min(cp);
% cpnorm = (cp + abs(min(cp)))./(4*max(abs(cp)));

% normalize Cp for otf
% cponorm = (cpo + abs(min(cpo)))./(4*max(abs(cpo)));


% calc d*
% flag(1:nx)=0;
kmax(1:nx) = nz/2;
umax(1:nx) = u0;
for i=1:nx
%     flag(i) = 0;
    for k=1:nzp-20
        if (k==1)
            dudz(i,k) = (u(i,k)-u(i,k+1))/(zpts(k)-zpts(k+1));
        else
            dudz(i,k) = (u(i,k-1)-u(i,k+1))/(zpts(k-1)-zpts(k+1));
        end
        
        if ( dudz(i,k) < 10^-6 && dudz(i,k)>0 && u(i,k) > 0.9*u0 )
            umax(i)=u(i,k);
            kmax(i) = k;
%             flag(i)=1;
        end
    end
end

kmax(1:nx) = 10;
dstar(1:nx)=0;
for i=1:nx
    sum=0.;
%     u0 = umax(i);
    u0 = u(i,kmax(i));
    for k=kmax(i):nz
        sum=sum-0.5*(1-u(i,k)/u0+1-u(i,k+1)/u0)*(zpts(k+1)-zpts(k));
        k=k-1;
    end
    dstar(i)=sum;
end

u0=1.;
% for i=1:1000
%     xi(i) = 9.6*10^5*xnu/(1000-1)*(i-1) + 10^5*xnu;
%     dblas(i) = 1.72*sqrt(xi(i)*xnu/u0); 
% end
% xi = xi/zpts(1);
dblas = 1.72*sqrt(xpts*xnu/u0);
dblas = dblas/zpts(1);


% % plot results

% scale x coordinates
x_scaled = xpts/zpts(1);
% scale vertical coordinates
z_scaled = zpts/zpts(1);
% scale spanwise coordinates
y_scaled = ypts/zpts(1);

% scale bl thicknesses
dstar = dstar/zpts(1);
% dstaro = dstaro/zpts(1);
% dstarb = dstarb/zpts(1);

% Cf
figure(1)
% plot(x_scaled,cf);
plot(x_sp,cf_sp,x_scaled,cf)
% plot(x_sp,cf_sp,x_scaled,cf,x_scaled,cfo,x_scaled,cfbl);
axis([x_scaled(1) x_scaled(nx) -0.007 0.01]);
xlabel('x/H');
ylabel('C_f');
% legend('Suction')
title('C_f plot');
legend('Spalart','UDNS');
% legend('Spalart','UDNS','UDNS OTF','Blasius IC');
saveas(1,'cf_otf_vs_qavg');

% Cp
figure(2)
% plot(x_scaled,cp);
plot(x_sp,cp_sp,x_scaled,cp);
% plot(x_sp,cp_sp,x_scaled,cp,x_scaled,cpo);
axis([x_scaled(1) x_scaled(nx) -0.01 0.6])
% axis([0.5 7.5 0. 0.5]);
xlabel('x/H');
ylabel('C_p');
legend('Spalart','UDNS')
% legend('Suction')
% legend('Spalart','UDNS','UDNS OTF');
saveas(2,'cp_otf_vs_qavg');

% contour plot of streamwise velocity
for k=1:nzp
    xz(:,k) = x_scaled;
end
for i=1:nx
    zx(i,:) = z_scaled;
end
for j=1:ny
    xy(:,j) = x_scaled;
end
for i=1:nx
    yx(i,:) = y_scaled;
end
figure(3)
contour(xz,zx,u,20);
colorbar;
xlabel('x/H');
ylabel('z/H');
title('time-avg contour plot u(x,z)');
saveas(3,'countour_u_otf_vs_qavg');

% dstar
figure(4)
% plot(x_scaled,dstar,xi,dblas)
plot(x_sp,dstar_sp,x_scaled,dstar,x_scaled,dblas)
% plot(x_sp,dstar_sp,x_scaled,dstar,x_scaled,dstarb,x_scaled,dstaro,xi,dblas);
title('Displacement thickness vs x');
xlabel('x/Y');
ylabel('\delta^*/Y');
legend('Spalart','UDNS','Blasius');
% legend('Spalart','UDNS','IC','UDNS OTF','Blasius');
saveas(4,'dstar_otf_vs_qavg');

% plot p(x,z->0)
figure(10)
plot(x_scaled,p(:,nzp),x_scaled,p(:,nz),x_scaled,p(:,nz-1),x_scaled,p(:,1))
legend(num2str(zpts(nzp)),num2str(zpts(nz)),num2str(zpts(nz-1)),num2str(zpts(1)));
ylabel('p(x)');
xlabel('x');
title('p(x) at diff z loc');
saveas(10,'p(x)_vs_z');

% plot vertical profiles of p(z)
figure(11)
% plot(p(170:190,:),z_scaled);
plot(p(1:nx/20:nx,:),z_scaled);
% plot(p(1:50:nx,nz-10:nzp),z_scaled(nz-10:nzp));
% plot(p(1:50:nx,1:nz-10),z_scaled(1:nz-10),'x');
xlabel('p(z)');
ylabel('z/H');
title('p(z) at diff x loc');
saveas(11,'p(z)_vs_x');

% u(z)
figure(5)
% plot(u(170:190,:),z_scaled);
plot(u(1:nx/20:nx,1:nzp),z_scaled(1:nzp));
xlabel('u(z) [cm/s]');
ylabel('z/H');
title('u(z) at diff x loc');
saveas(5,'u(z)_vs_x');

% u(x)
figure(6)
plot(x_scaled,u(:,5:1:nz));
ylabel('u(x) [cm/s]');
xlabel('x/H');
title('u(x) at diff z loc');
saveas(6,'u(x)_vs_z');

% w(x)
figure(7)
plot(x_scaled,w(:,1:4:nz));
ylabel('w(x) [cm/s]');
xlabel('x/H');
title('w(x) at diff z loc');
saveas(7,'w(x)_vs_z');

% w(z)
figure(8)
% plot(w(170:190,:),z_scaled);
plot(w(1:nx/20:nx,1:nzp),z_scaled(1:nzp));
xlabel('w(z) [cm/s]');
ylabel('z/H');
title('w(z) at diff x loc');
saveas(8,'w(z)_vs_x');

% v(z)
figure(88)
% plot(w(170:190,:),z_scaled);
plot(v(1:nx/20:nx,1:nzp),z_scaled(1:nzp));
xlabel('v(z) [cm/s]');
ylabel('z/H');
title('v(z) at diff x loc');
saveas(88,'v(z)_vs_x');

% temp(z)
figure(55)
% plot(u(170:190,:),z_scaled);
plot(temp(1:nx/20:nx,1:nzp),z_scaled(1:nzp));
xlabel('temp(z) [cm/s]');
ylabel('z/H');
title('temp(z) at diff x loc');
saveas(55,'temp(z)_vs_x');


% w(x)
figure(77)
plot(x_scaled,temp(:,nz-10:nz));
ylabel('temp(x) [cm/s]');
xlabel('x/H');
title('temp(x) at diff z loc');
saveas(77,'temp(x)_vs_z');


% spanwise plots
wy(:,:) = win(:,:,nzp-20);
uy(:,:) = uin(:,:,nzp-20);
vy(:,:) = vin(:,:,nzp-20);
% 
% 
figure(30)
% contour(xy,yx,wy);
plot(y_scaled,wy(nx/3:nx/10:nx,:))
% colorbar;
xlabel('y/H');
ylabel('w(y)');
title('w(y) at different x locations');
% ylabel('y/H');
% title('time-avg contour plot w(x,y,nzp-25)');
saveas(30,'w(y)_vs_x');
% 
figure(31)
% contour(xy,yx,wy);
plot(y_scaled,vy(nx/3:nx/10:nx,:))
% colorbar;
xlabel('y/H');
ylabel('v(y)');
title('v(y) at different x locations');
% ylabel('y/H');
% title('time-avg contour plot w(x,y,nzp-25)');
saveas(31,'v(y)_vs_x');

figure(32)
% contour(xy,yx,wy);
plot(y_scaled,uy(nx/3:nx/10:nx,:))
% colorbar;
xlabel('y/H');
ylabel('u(y)');
title('u(y) at different x locations');
% ylabel('y/H');
% title('time-avg contour plot w(x,y,nzp-25)');
saveas(32,'u(y)_vs_x');


% boundary conditions
figure(101)
plot(x_scaled,u(:,1)-1.,x_scaled,u(:,nzp),x_scaled,w(:,1),x_scaled,w(:,nzp));
xlabel('x/Y');
ylabel('velocity [cm/s]');
title('Boundary conditions');
legend('u(Y) - ubar(Y)','u(y=0)','w(Y)','w(y=0)');

% plot dzplus to get a sense of resolution
% figure(102)
% plot(x_scaled,dzplus);
% xlabel('x');
% ylabel('dz^+')