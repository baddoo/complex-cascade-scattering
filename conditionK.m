omega5 = .5; he = .5; w = 1; sigma5 = 3*pi/4; d =1;
mu = [1 0 0]; s = sqrt(d^2 + he^2);


gamma = 1*(linspace(-5,5,5e2) + 1i*linspace(-5,5,5e2)');

AAData.omega5 = omega5;
AAData.w = w;
ADData.spac(1) = he;
ADData.spac(2) = d;
AAData.sigma5 = sigma5;
ADData.mu = mu;

zeta = mysqrt(omega5*w,gamma);
J = Jfun(gamma,ADData,AAData);
approx = zeta.*tan(he*zeta)./...
        (4*pi*(1-cos(d*gamma + sigma5)./cos(he*zeta)));
approx = zeta.*sin(he*zeta)./...
        (4*pi*(cos(he*zeta)-cos(d*gamma + sigma5)));  
approx = zeta/1i.*(exp(1i*he*zeta)-exp(-1i*he*zeta))./...
        (4*pi*(exp(1i*he*zeta)+exp(-1i*he*zeta)-exp(1i*(d*gamma + sigma5))-exp(-1i*(d*gamma+sigma5))));     
approx = zeta/1i.*(exp(2i*he*zeta)-1)./...
        (4*pi*(exp(2i*he*zeta)+1-exp(1i*(he*zeta + d*gamma + sigma5))-exp(-1i*(-he*zeta + d*gamma+sigma5))));     
%approx = (exp(1i*(d*gamma+sigma)) + exp(-1i*(d*gamma + sigma)))...
%       ./(exp(1i*he*zeta)  + exp(-1i*he*zeta)); 
%approx = (exp(1i*(2*d*gamma - he*zeta +sigma)) + exp(-1i*(he*zeta + sigma)))...
%        ./(exp(1i*d*gamma)  + exp(-1i*(2*he*zeta-d*gamma))); 
%approx = (exp(1i*(d*gamma - he*zeta +sigma)) + exp(-1i*(d*gamma + he*zeta + sigma)))...
%        ./(1  + exp(i*2*he*zeta));     
%approx = (./(1+exp(1i*2*zeta));    
subplot(1,2,1)
pcolor(real(gamma),imag(gamma),(abs(J))); shading interp;
caxis([0,10])
subplot(1,2,2)
pcolor(real(gamma),imag(gamma),(abs(approx))); shading interp;
caxis([0,10])
colormap jet

% Condition argument
zeta = @(x) mysqrt(omega5*w,x);
sigma = sigma5;
den = @(x) cos(he*zeta(x)) - cos(d*x + sigma);    
dden= @(x) he*x./zeta(x).*sin(he*zeta(x)) + d*sin(d*x + sigma);
num = @(x) zeta(x).*sin(he*zeta(x)) + (mu(1) + x*mu(2) + x.^2*mu(3)).*den(x);
dnum= @(x) -x./zeta(x).*sin(he*zeta(x)) - he*x.*cos(he.*zeta(x)) + (mu(1) + x*mu(2) + x.^2*mu(3)).*dden(x) + (mu(2) + 2*x.*mu(3)).*den(x);
argFun = @(x) dnum(x)./num(x);
%apDNum = @(x) ;
p = @(x) mu(1) -1i* x*mu(2) - x.^2*mu(3);
Dp= @(x) -1i*mu(2) - 2*x*mu(3);
denR = @(x) (exp(2i*he*zeta(x))+1 - exp(1i*(he*zeta(x)+ d*x + sigma))-exp(-1i*(-he*zeta(x)+ d*x + sigma)));
%denR = @(x)  exp(1i*(2*he*zeta(x)+d*x))+exp(1i*d*x) - exp(1i*(he*zeta(x)+ 2*d*x + sigma))-exp(-1i*(-he*zeta(x)+ sigma));
denRP = @(x) (exp(1i*(2*he*zeta(x)+d*x))+exp(1i*d*x) - exp(1i*(he*zeta(x)+ 2*d*x + sigma))-exp(-1i*(-he*zeta(x)+ sigma)));
denRP2 = @(x) (exp(1i*(he*zeta(x)+d*x))+exp(1i*(d*x-he*zeta(x))) - exp(1i*(2*d*x + sigma))-exp(-1i*( sigma)));

Num1R = @(x) zeta(x)/1i.*(exp(2i*he*zeta(x))-1);
Num1RP = @(x) zeta(x)/1i.*(exp(1i*(2*he*zeta(x)+d*x))-exp(1i*d*x));
Num1RP2= @(x) zeta(x)/1i.*(exp(1i*(he*zeta(x)+d*x))-exp(1i*(d*x-he*zeta(x))));

dNum1R  = @(x) -x./zeta(x)/1i.*(exp(2i*he*zeta(x))-1) - he*x.*(exp(2i*he*zeta(x))+1);
dNum1RP = @(x) -x./zeta(x)/1i.*(exp(1i*(2*he*zeta(x)+d*x))-exp(1i*d*x)) - he*x.*(exp(1i*(2*he*zeta(x)+d*x))+exp(1i*d*x));
dNum1RP2= @(x) -x./zeta(x)/1i.*(exp(1i*(he*zeta(x)+d*x))-exp(1i*(d*x-he*zeta(x)))) - he*x.*(exp(1i*(he*zeta(x)+d*x))+exp(1i*(d*x-he*zeta(x))));

dDenR  = @(x)  (he*x./zeta(x)/1i.*(exp(2i*he*zeta(x))-1) + d/1i*(exp(1i*(he*zeta(x)+d*x + sigma))-exp(-1i*(-he*zeta(x)+d*x + sigma))));
dDenRP = @(x) he*x./zeta(x)/1i.*(exp(1i*(2*he*zeta(x)+d*x))-exp(1i*d*x)) + d/1i*(exp(1i*(he*zeta(x)+2*d*x + sigma))-exp(-1i*(-he*zeta(x) + sigma)));
dDenRP2= @(x) he*x./zeta(x)/1i.*(exp(1i*(he*zeta(x)+d*x))-exp(1i*(d*x-he*zeta(x)))) + d/1i*(exp(1i*(2*d*x + sigma))-exp(-1i*(+ sigma)));

apNum  = @(x) Num1RP2(x) + p(x).*denRP2(x);
apDNum = @(x) dNum1RP2(x) + p(x).*dDenRP2(x) + Dp(x).*denRP2(x);
approx = @(x) apDNum(x)./apNum(x);
%approx = @(x) (x.*(1+zeta(x)*he.*cot(he*zeta(x))))./-zeta(x).^2;

subplot(1,2,1)
PhasePlot(gamma,KNumlogD(gamma,ADData,AAData)); axis on
%pcolor(real(gamma),imag(gamma),(imag(argFun(gamma)))); shading interp;
%caxis([0,1e2])
subplot(1,2,2)
pcolor(real(gamma),imag(gamma),imag(KNumlogD(gamma,ADData,AAData))); shading interp;
%caxis([0,1e2])
colorbar
norm(abs(argFun(gamma)-KlogD(gamma,ADData,AAData)),'inf')