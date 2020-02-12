addpath('cfe')

x = linspace(-20,20,1000);
y = linspace(-30,30,1000).';
z = x + 1i*y;

d = .4; s = .5;
k=3; w=1;
mu=.03; sigma =2;

f = @(x) (mysqrt(k*w,x).*sin(s*mysqrt(k*w,x))+ mu*(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma))) ...
          ./(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma));

      r = linspace(0,1,100); th = linspace(-pi,pi,100);
      
  Z = r.*exp(1i*th.');
  a = .5;
  fl = @(zVar) 1/2i/pi*log((zVar-a)./abs(a)./(zVar - 1./conj(a)));
PhasePlot(z,f(z),'d');% shading interp
hold on
%plot3([0,10],[0,10],[0,0],'LineWidth',3)
plot([0,-10],[0,10],'LineWidth',3)
hold off
%colorbar