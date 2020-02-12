addpath('cfe')
% CHANGE BACK TO WAYPOINTS
d = .5; s = 1;
k=20*(1+0.00i); w=1;
mu=1; sigma =3;
zeta = @(xVar) mysqrt(k*w,xVar);
f = @(x) (mysqrt(k*w,x).*sin(s*mysqrt(k*w,x))+ mu*(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma))) ...
          ./(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma));
fn  = @(x) (zeta(x).*sin(s*zeta(x))+ mu*(cos(s*zeta(x)) - cos(d*x + sigma)));
fd = @(x) (cos(s*mysqrt(k*w,x)) - cos(d*x + sigma));
fnp = @(x) -x./zeta(x).*sin(s*zeta(x))-s*x.*cos(s*zeta(x)) + mu*(s*x./zeta(x).*sin(s*zeta(x)) + d*sin(d*x+sigma));
argFun = @(xVar) fnp(xVar)./fn(xVar);

%r = 5;
l=2*k/d;h=2;
%argFunThetInt = @(thetVar) 1i*r.*exp(1i*thetVar).*argFun(r.*exp(1i*thetVar));
ways = [(l-1i*h)/2,(l+1i*h)/2,(-l+1i*h)/2];
allWays = [-(l+1i*h)/2,ways,-(l+1i*h)/2];
Npr = (real(integral(argFun,-(l+1i*h)/2,-(l+1i*h)/2,'Waypoints',ways,'ArrayValued',true)/(2i*pi)));
N = round(Npr);
N3 = permute(1:N,[1,3,2]);
polyFun = @(xVar) bsxfun(@power,xVar,N3);
newFun = @(zVar) polyFun(zVar).*argFun(zVar);
tic
p = integral(newFun,-(l+1i*h)/2,-(l+1i*h)/2,'Waypoints',ways,'ArrayValued',true)/(2i*pi);

e = ones(1,N+1);

for l = 1 :N
    j = (1:l).';
    e(l+1) = 1./l.*sum(permute((-1).^(j-1).'.*e(l-j+1),[1,3,2]).*p(:,:,(1:l)),3);
end

myroots = roots((-1).^(1:N+1).*e);
toc
nx = 500; ny = 500;
xp = linspace(-15,15,nx);
yp = linspace(-30,30,ny);
[X,Y]=meshgrid(xp,yp);
n = 0:200;
Rsols1  = mu - 1i/(d-1i*s)*lambertw(+n,-1i*mu*(d-1i*s)*exp(-1i*(mu*(d-1i*s)+sigma)));
Rsols2 = -mu - 1i/(d+1i*s)*lambertw(-n,+1i*mu*(d+1i*s)*exp(+1i*(mu*(d+1i*s)-sigma)));

Rsols1a = mu + (-(1+4*(0:300))*pi+2*(mu+sigma))./(2*(d-1i*s));
a  = (d+1i*s)*2*pi/(d^2+s^2);
b  = 1i*(d+1i*s)/(d^2+s^2);
c  = (pi-sigma+1i*log(2*pi/mu/sqrt(d^2+s^2))-atan(s/d))*(d+1i*s)./(d^2+s^2);

a1 = (-d+1i*s)*2*pi/(d^2+s^2);
b1 =  -1i*(-d+1i*s)/(d^2+s^2);
c1 = (-pi+sigma-1i*log(2*pi/mu/sqrt(d^2+s^2))-atan(s/d))*(-d+1i*s)./(d^2+s^2);

n1 = n;
Rsols1b = a*n+b*log(n)+c;
Rsols2b = a1*n+b1*log(n)+c1;
Rsols3b = conj(a1)*n+conj(b1)*log(n)+conj(c1);
Rsols4b = conj(a)*n+conj(b)*log(n)+conj(c);

figure(1)
semilogy(abs(Rsols2./Rsols2b-1));
figure(2)
plot(abs(Rsols2-Rsols2b))
hold on
%plot(abs(Rsols1b))
hold off

%Rsols1(abs(Rsols1)<r/2)=[];
%Rsols3(abs(Rsols3)<r/2)=[];
Rsols4 = conj(Rsols1);
Rsols3 = conj(Rsols2);

% figure(2)
% semilogy(abs(f((Rsols1))))

%return
figure(1)
Z = X + 1i*Y;
PhasePlot(Z,fn(Z)/10,'d');
%return
hold on
plot(myroots,'x')

%plot(r.*exp(1i*linspace(-pi,pi)),'k','LineWidth',3)
axis on

tic 
sols1a=zeros(2,numel(Rsols1));
sols2a=zeros(2,numel(Rsols2));
sols3a=zeros(2,numel(Rsols3));
sols4a=zeros(2,numel(Rsols4));
sols1aM = zeros(1,numel(Rsols1));

for nL = 1:numel(Rsols1)
 sols1a(:,nL) = findRoots(s,d,sigma,mu,k,w,[real(Rsols1(nL));imag(Rsols1(nL))]);
 sols2a(:,nL) = findRoots(s,d,sigma,mu,k,w,[real(Rsols2(nL));imag(Rsols2(nL))]);
 sols3a(:,nL) = findRoots(s,d,sigma,mu,k,w,[real(Rsols3(nL));imag(Rsols3(nL))]);
 sols4a(:,nL) = findRoots(s,d,sigma,mu,k,w,[real(Rsols4(nL));imag(Rsols4(nL))]);
end
toc

tic 
sols1ac=zeros(1,numel(Rsols1));
sols2ac=zeros(1,numel(Rsols2));
sols3ac=zeros(1,numel(Rsols3));
sols4ac=zeros(1,numel(Rsols4));

for nL = 1:numel(Rsols1)
 sols1ac(nL) = findRoots2(s,d,sigma,mu,k,w,Rsols1(nL));
 sols2ac(nL) = findRoots2(s,d,sigma,mu,k,w,Rsols2(nL));
 sols3ac(nL) = findRoots2(s,d,sigma,mu,k,w,Rsols3(nL));
 sols4ac(nL) = findRoots2(s,d,sigma,mu,k,w,Rsols4(nL));
end
toc

tic 

sols1ab = findRoots2(s,d,sigma,mu,k,w,Rsols1);
sols2ab = findRoots2(s,d,sigma,mu,k,w,Rsols2);
sols3ab = conj(sols2ab);
sols4ab = conj(sols1ab);

toc


%max(max(abs(sols1-sols1a)))
%max(abs(f(myroots)))
%max(abs(f(complex(sols1(1,:),sols1(2,:)))))
%max(max(abs(f(complex(sols1a(1,:),sols1a(2,:))))))

% %return
% sols2 = findRoots(s,d,sigma,mu,k,w,[real(Rsols2);imag(Rsols2)]);
% sols3 = findRoots(s,d,sigma,mu,k,w,[real(Rsols3);imag(Rsols3)]);
% sols4 = findRoots(s,d,sigma,mu,k,w,[real(Rsols4);imag(Rsols4)]);

comSols1 = complex(sols1a(1,:),sols1a(2,:));
comSols2 = complex(sols2a(1,:),sols2a(2,:));
comSols3 = complex(sols3a(1,:),sols3a(2,:));
comSols4 = complex(sols4a(1,:),sols4a(2,:));

plot([comSols1;Rsols1],'--')
plot(comSols1,'x')
plot(Rsols1,'o')

plot([comSols2;Rsols2],'--')
plot(comSols2,'x')
plot(Rsols2,'o')

plot([comSols3;Rsols3],'--')
plot(comSols3,'x')
plot(Rsols3,'o')

plot([comSols4;Rsols4],'--')
plot(comSols4,'x')
plot(Rsols4,'o')

plot(allWays,'k','LineWidth',3)
axis equal
axis([xp(1),xp(end),yp(1),yp(end)])
hold off
%%
figure(2) 
loglog(n,abs(comSols1-Rsols1b),'o')
hold on
loglog(n,abs(comSols2-Rsols2b),'s')
hold off
axis tight
xlabel('$n$')
return
%% Plot solutions
phi = linspace(0,pi*2,100).'; phi(1) = []; phi(end) = [];
phi= pi/3;
R = exp(linspace(1,15));

bigZ = R.*exp(1i*phi);
zeta = mysqrt(k*w,bigZ);
test1 = sin(s*zeta)./(-exp(s*bigZ)/2i);
%test1 = cos(s*zeta)./(exp(s*bigZ)/2);
%test1 = cos(d*bigZ+sigma)./(exp(-1i*d*bigZ-1i*sigma)/2);

figure(2)
%pcolor(real(bigZ),imag(bigZ),real(zeta)); shading interp
%return
%plot(abs(zeta-1i*sqrt(R).*exp(1i*phi.'/2)))
%plot(R,imag(test1));% -  1/2i*exp(-1i*h*bigZ)))
%hold on
plot(R,real(test1));%-1/2i*exp(h*R*cos(phi)+1i*h*R*sin(phi))))
hold on
plot(R,imag(test1))
%plot(abs(sqrt(R).*exp(1i*phi.'/2)))
%plot(abs(sqrt(bigZ)).')
%plot(abs(a).')
hold off


