function Au = Au(X,Y,Gamma,Zeta,spac,Sigma)

s = spac(1); d = spac(2); % Spacing
Kd = 4*pi*(cos(s*Zeta) - cos(d*Gamma + Sigma));
Au = -2*pi*exp(1i*(d*Gamma + Sigma)).*cos(Zeta.*Y).*exp(-1i*Gamma.*X)./Kd;

end