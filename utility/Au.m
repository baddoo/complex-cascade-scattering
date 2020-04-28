function Au = Au(X,Y,Gamma,Zeta,ADData)

s = ADData.spac(1); d = ADData.spac(2); % Spacing
sigma=data.sigma; % Inter-blade phase angle
Kd = 4*pi*(cos(s*Zeta) - cos(d*Gamma + sigma));
Au = -2*pi*exp(1i*(d*Gamma + sigma)).*cos(Zeta.*Y).*exp(-1i*Gamma.*X)./Kd;

end