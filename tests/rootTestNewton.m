k = 20; s = .4; w = .3; sigma = 2; d = 1/2;
mu = 1;

zeta = @(xVar) mysqrt(k*w,xVar);

%f = @(x) (mysqrt(k*w,x).*sin(s*mysqrt(k*w,x))+ mu*(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma))) ...
%          ./(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma));
den = @(x) cos(s*zeta(x)) - cos(d*x + sigma);    
dden= @(x) s*x./zeta(x).*sin(s*zeta(x)) + d*sin(d*x + sigma);
num = @(x) zeta(x).*sin(s*zeta(x));
dnum= @(x) -x./zeta(x).*sin(s*zeta(x)) - s*x.*cos(s.*zeta(x));
f = @(x) num(x)./den(x) + mu*x;
df= @(x) (dnum(x).*den(x) - num(x).*dden(x))./den(x).^2 + mu;
fn  = @(x) (zeta(x).*sin(s*zeta(x))+ mu*(cos(s*zeta(x)) - cos(d*x + sigma)));
%fd = @(x) (cos(s*mysqrt(k*w,x)) - cos(d*x + sigma));
fnp = @(x) -x./zeta(x).*sin(s*zeta(x))-s*x.*cos(s*zeta(x)) + mu*(s*x./zeta(x).*sin(s*zeta(x)) + d*sin(d*x+sigma));


kV = permute(0:5,[1,3,2]);
argFun = @(xVar) dnum(xVar)./num(xVar).*xVar.^kV;
argFun = @(xVar) df(xVar)./f(xVar).*xVar.^kV;
%argFun = @(xVar) 1./xVar;

R = pi*permute(1:5,[1,4,3,2]);
intFun = @(th) 1i*R.*exp(1i*th).*argFun(R.*exp(1i*th));
tic
I13 = integral(intFun,0,2*pi,'ArrayValued',true)/(2i*pi);
I1 = permute(I13,[3,4,1,2])
toc
N = 1e3;
tInt = linspace(0,2*pi,N+1); tInt(end) = [];
tic
I23 = 2*pi/N*sum(intFun(tInt),2)/(2i*pi);
I2 = permute(I23,[3,4,1,2])
toc
err = norm(I1-I2,'inf')
return
%argFun = @(xVar) dnum(xVar)./num(xVar);
%f = @(x) 1./x.^3;
%fp = @(x) -3./x.^4;
R = 50;

th = linspace(0,2*pi);
xp = linspace(-R,R,1e4);
xu = xp + 3i;
xd = xp - 3i;
%plot(th,abs(argFun(R*exp(1i*th))));
pF = argFun(xu) - argFun(xd);
plot(xp,abs(pF))
%semilogy(xp,abs(argFun(xp).*xp-1));
return

R = 500;
integral(@(th) 1i*(R*exp(1i*th)).*argFun(R*exp(1i*th)),0,2*pi)/(2i*pi)
tInt = linspace(0,2*pi,3000);
trapz(tInt, 1i*(R*exp(1i*tInt)).*argFun(R*exp(1i*tInt)))/(2i*pi)
return
tol = 1e-4;
x0 = 100*rand(1);
err = inf; 
j = 0; % counter
tic
xn = x0;
while err>tol
     fn=f(xn); %Calculating the value of function at x0
     fn_der=fp(xn); %Calculating the value of function derivative at x0
  xn=xn-fn/fn_der; % The Formula
err = abs(f(xn));
j = j+1;
end
disp(['The time for the Newton solver was ',num2str(toc),' seconds.']);
tic
options = optimset('Display','off');
tic
x0F = fsolve(f,x0,options);
disp(['The time for fsolve was ',num2str(toc),' seconds.']);
[abs(f(xn)),abs(f(x0F))]
