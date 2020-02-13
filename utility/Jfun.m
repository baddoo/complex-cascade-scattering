function J = Jfun(gamma,ADData,AAData)

omega5 = AAData.omega5;
w = AAData.w;
he = ADData.spac(1);
d = ADData.spac(2);
sigma5 = AAData.sigma5;
mu = ADData.mu;

%J = zeros(size(gamma./omega5));

J= mysqrt(omega5*w,gamma).*sin(he*mysqrt(omega5*w,gamma))./(4*pi*(cos(he*mysqrt(omega5*w,gamma)) - cos(d*gamma + sigma5)));
        
end