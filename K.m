function K = K(gamma,ADData,AAData)

omega5 = AAData.omega5;
w = AAData.w;
he = ADData.spac(1);
d = ADData.spac(2);
sigma5 = AAData.sigma5;
mu = ADData.mu;

K= (mysqrt(omega5*w,gamma).*sin(he*mysqrt(omega5*w,gamma)) + (mu(1)-1i*gamma*mu(2)-mu(3)*gamma.^2).*(cos(he*mysqrt(omega5*w,gamma)) - cos(d*gamma + sigma5))) ...
            ./(4*pi*(cos(he*mysqrt(omega5*w,gamma)) - cos(d*gamma + sigma5)));
%K= (mysqrt(omega5*w,gamma).*sin(he*mysqrt(omega5*w,gamma))./(cos(he*mysqrt(omega5*w,gamma)) - cos(d*gamma + sigma5))/(4*pi) + (mu(1)-1i*gamma*mu(2)-mu(3)*gamma.^2)/4/pi);

end