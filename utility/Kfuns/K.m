function K = K(gamma,ADData,AAData)

omega = AAData.omega;
w = AAData.w;
s = ADData.spac(1);
d = ADData.spac(2);
sigma = AAData.sigma;
mu = ADData.mu;

K= (mysqrt(omega*w,gamma).*sin(s*mysqrt(omega*w,gamma)) + (mu(1)-1i*gamma*mu(2)-mu(3)*gamma.^2).*(cos(s*mysqrt(omega*w,gamma)) - cos(d*gamma + sigma))) ...
            ./(4*pi*(cos(s*mysqrt(omega*w,gamma)) - cos(d*gamma + sigma)));

end