function K = regK(gamma,ADData,AAData)

% Extract data from structures
omega = AAData.omega;
w = AAData.w;
s = ADData.spac(1);
d = ADData.spac(2);
Sigma = AAData.Sigma;
mu = ADData.mu; 

zeta = @(x) mysqrt(omega*w,x);
zg = zeta(gamma);

regJ = zg/1i.*(exp(2i*s*zg)-1)./...
        (4*pi*(exp(2i*s*zg)+1-exp(1i*(s*zg + d*gamma + Sigma))-exp(-1i*(-s*zg + d*gamma+Sigma))));     

K = regJ + (mu(1) - 1i*gamma*mu(2) - gamma.^2*mu(3))/(4*pi);    

end