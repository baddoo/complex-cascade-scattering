function Ad = Ad(X,Y,Gamma,Zeta,spac,Sigma,PM,M,type)
%AD is a function used to evaluate the solution.

s = spac(1); d = spac(2); % Spacing
Kd = 4*pi*(cos(s*Zeta) - cos(d*Gamma + Sigma));

if strcmp(type,'vvelocity')
Ad = -2*pi*Zeta.*sin(Zeta.*(Y - s)).*exp(-1i*Gamma.*X)./Kd.*exp(1i*PM*M^2*X);    
else
    
Ad = 2*pi*cos(Zeta.*(Y - s)).*exp(-1i*Gamma.*X)./Kd;

    if strcmp(type,'pressure')
        Ad = -1i*(PM - Gamma).*exp(1i*PM*M^2*X).*Ad;
    elseif strcmp(type,'hvelocity')
        Ad  = -1i*Gamma.*exp(1i*PM*M^2*X).*Ad;
    end
    
end



end