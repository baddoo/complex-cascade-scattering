function K = regK(gamma,ADData,AAData)

% Extract data from structures
omega5 = AAData.omega5;
w = AAData.w;
he = ADData.spac(1);
d = ADData.spac(2);
sigma = AAData.sigma5;
mu = ADData.mu; 

zeta = @(x) mysqrt(omega5*w,x);
zg = zeta(gamma);

regJ = zg/1i.*(exp(2i*he*zg)-1)./...
        (4*pi*(exp(2i*he*zg)+1-exp(1i*(he*zg + d*gamma + sigma))-exp(-1i*(-he*zg + d*gamma+sigma))));     

K = regJ + (mu(1) - 1i*gamma*mu(2) - gamma.^2*mu(3))/(4*pi);    

end