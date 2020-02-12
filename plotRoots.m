addpath('cfe')
% CHANGE BACK TO WAYPOINTS

d = 0.5; s = 1;
k=2*(1+0.00i); w=1;
mu=5; sigma =3;
nRoots = 250;

F= @(x) (mysqrt(k*w,x).*sin(s*mysqrt(k*w,x)) + mu*(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma))) ...
            ./(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma));

n = 1:floor(nRoots/nRoots):nRoots;
Rsols1  = mu - 1i/(d-1i*s)*lambertw(+n,-1i*mu*(d-1i*s)*exp(-1i*(mu*(d-1i*s)+sigma)));
Rsols2 = -mu - 1i/(d+1i*s)*lambertw(-n,+1i*mu*(d+1i*s)*exp(+1i*(mu*(d+1i*s)-sigma)));

a  = (d+1i*s)*2*pi/(d^2+s^2);
b  = 1i*(d+1i*s)/(d^2+s^2);
c  = (pi-sigma+1i*log(2*pi/mu/sqrt(d^2+s^2))-atan(s/d))*(d+1i*s)./(d^2+s^2);

a1 = (-d+1i*s)*2*pi/(d^2+s^2);
b1 =  -1i*(-d+1i*s)/(d^2+s^2);
c1 = (-pi+sigma-1i*log(2*pi/mu/sqrt(d^2+s^2))-atan(s/d))*(-d+1i*s)./(d^2+s^2);

Rsols1b = a*n+b*log(n)+c;
Rsols2b = a1*n+b1*log(n)+c1;

sols1ab = findRoots2(s,d,sigma,mu,k,w,Rsols1b);
sols2ab = findRoots2(s,d,sigma,mu,k,w,Rsols2b);

figure(1) 
semilogy(n,abs(sols1ab-Rsols1b),'o')
hold on
semilogy(n,abs(sols2ab-Rsols2b),'s')
hold off
axis tight
xlabel('$n$')


figure(2)

plot(