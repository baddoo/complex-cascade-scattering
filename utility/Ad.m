function Ad = Ad(X,Y,Gamma,Zeta,spac,Sigma)

s = spac(1); d = spac(2); % Spacing
Kd = 4*pi*(cos(s*Zeta) - cos(d*Gamma + Sigma));
Ad = 2*pi*cos(Zeta.*(Y - s)).*exp(-1i*Gamma.*X)./Kd;

end