function Ad = Ad(X,Y,Gamma,Zeta,ADData)

s = ADData.spac(1); d = ADData.spac(2); % Spacing
sigma=data.sigma; % Inter-blade phase angle
Kd = 4*pi*(cos(s*Zeta) - cos(d*Gamma + sigma));
Ad = 2*pi*cos(Zeta.*(Y - s)).*exp(-1i*Gamma.*X)./Kd;

end