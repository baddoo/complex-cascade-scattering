function Au = Au(X,Y,Gamma,Zeta,spac,Sigma,PM,M,type)

s = spac(1); d = spac(2); % Spacing
Kd = 4*pi*(cos(s*Zeta) - cos(d*Gamma + Sigma));

if strcmp(type,'vvelocity')
Au = 2*pi*exp(1i*(d*Gamma + Sigma)).*Zeta.*sin(Zeta.*Y).*exp(-1i*Gamma.*X)./Kd.*exp(1i*PM*M^2*X);    
else
    
Au = -2*pi*exp(1i*(d*Gamma + Sigma)).*cos(Zeta.*Y).*exp(-1i*Gamma.*X)./Kd;

    if strcmp(type,'pressure')
        Au = 1i*(PM - Gamma).*exp(1i*PM*M^2*X).*Au;
    elseif strcmp(type,'hvelocity')
        Au  = -1i*Gamma.*exp(1i*PM*M^2*X).*Au;
    end
    
end


end