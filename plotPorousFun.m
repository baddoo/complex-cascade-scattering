addpath('cfe')
% CHANGE BACK TO WAYPOINTS
d = 0.5; s = 1;
k=2*(1+0.00i); w=1;
mu=20; sigma =3;
zeta = @(xVar) mysqrt(k*w,xVar);
f = @(x) (mysqrt(k*w,x).*sin(s*mysqrt(k*w,x))+ mu*(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma))) ...
          ./(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma));
fn  = @(x) (zeta(x).*sin(s*zeta(x))+ mu*(cos(s*zeta(x)) - cos(d*x + sigma)));
fd = @(x) (cos(s*mysqrt(k*w,x)) - cos(d*x + sigma));
fnp = @(x) -x./zeta(x).*sin(s*zeta(x))-s*x.*cos(s*zeta(x)) + mu*(s*x./zeta(x).*sin(s*zeta(x)) + d*sin(d*x+sigma));

figure(1)
nx = 400;
ny = 200;
nd = 20; ns = 20;
x = nd*d*linspace(-1,1,nx);
y = ns*s*linspace(-1,1,ny);
z = x+ 1i*y.';
P = PhasePlot(z,fn(z),'d'); shading interp
box on
axis on
%return
ax = gca;
%view([0,-90])
set(ax,'XAxisLocation','origin')
set(ax,'YAxisLocation','origin')
colorbar