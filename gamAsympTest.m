
%% Define parameter
d = 0.1; s = 1;
k=3*(1+0.00i); w=1;
mu=.1; sigma =1;

n = 1:1e6; % Number of terms we take the product over

nz = 45; % Number of terms in line
z = (-1-1i)*linspace(1000,5e3,nz).'; % Defines a line tending to infinity

a1 =  (d+1i*s)/(d^2+s^2)*2*pi; 
a2 =  (-d+1i*s)/(d^2+s^2)*2*pi;
a3 = a1;
a4 = a2;
b1 = 1i*(d+1i*s)/(d^2+s^2);
b2 = -b1*a2/a1;
c1  = ( pi-sigma+1i*log(2*pi/mu/sqrt(d^2+s^2))-atan(s/d))*( d+1i*s)./(d^2+s^2);
c2 =  (-pi+sigma-1i*log(2*pi/mu/sqrt(d^2+s^2))-atan(s/d))*(-d+1i*s)./(d^2+s^2);
c3 = -a1*sigma/2/pi;
c4 =  a2*sigma/2/pi;
d1 = 10;
d2 = 10;

fp = (sigma - 2*pi*n)/sqrt(d^2+s^2);
fm = (sigma + 2*pi*n)/sqrt(d^2+s^2);
chie = atan(d/s);
zeta = @(zVar) mysqrt(k*w,zVar);

lambdap = -fp*sin(chie)+cos(chie)*zeta(fp);
lambdam = -fm*sin(chie)+cos(chie)*zeta(fm); 

%newP = prod((1-z./(a1.*n + b1*log(n) + c1)).*(1-z./(a2*n + b2*log(n) + c2)) ...
%          ./(1-z./(a3.*n + c3))./(1-z./(a4*n + c4))    ,2);
%finProd2 = prod((1 - (b1*log(n))./(z-(a1*n+c1))).*(1 - (b2*log(n))./(z-(a2*n+c2))),2);
finProd = prod((1 - z./(a1*n+b1*log(n)+c1+d1./n)).*(1 - z./(a2*n+b2*log(n)+c2+d2./n))...
                ./((1-z./lambdap).*(1-z./lambdam)),2);
%finProd = prod(1 + b1*log(n)./(z-a1*n),2);
%finProd1 = sum(log(abs((1 - b1*log(n)./(z-(a1*n))).*(1 - 0*(b2*log(n))./(z-(a2*n+c2))))),2);
% 
% z1 = (1+1i)*linspace(-10,10,nz).';
% finProd = prod(1+1i./(n-z1),2);
% sum1 = cumsum(log(abs(1+1i./(n-z1))),2);


f = z.^(-c1./a1-c2./a2+c3./a3+c4./a4);

%%
figure(1)
loglog(abs(z),abs(finProd))

%%
% This figure checks the asymptotic behaviour of the second infinite
% product. In Nigel's paper it tends to 1. For us I suspect that it has
% some mixture of logarithmic and algebraic growth.

figure(2)

loglog(abs(z),((abs(finProd2))))

return

%% Verifies our description of the asymptotic behaviour of the gamma function
z1 = linspace(1,1e3,nz);

a = -4; b = 1;
myP = gamma((b-z1)./a);
aP = sqrt(2*pi).*(-z1/a).^(-z1/a+b/a-.5).*exp(z1/a);
  
%A = 1; B = 1;
%myP = gamma(A*z+B);
%aP = sqrt(2*pi).*exp(-A*z).*(A*z).^(A*z+B-.5);
%myP = prod((1-z./(a.*n + b)).*exp(z./(a.*n)),2);

%C = -b./sqrt(2*pi).*gamma(b./a).*((-a).^(b/a-.5));
%aP = C.*exp(z*(eulergamma-1-log(-a))/a).*z.^(z/a-b/a-.5);

figure(3)
semilogy(z1,abs(myP./aP-1))
% This plot should be decreasing.