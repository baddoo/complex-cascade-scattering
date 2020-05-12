function AdR = AdR(X,Y,Lambda,Zeta,SQRT,spac,PM,M,type)

s = spac(1); d = spac(2); % Spacing
del = sqrt(s^2 + d^2);
sgn = sign(imag(Lambda(end)));

AdRDen = 2*del*sin(s*Zeta).*SQRT;

if strcmp(type,'vvelocity')
AdRNum = -sgn*Zeta.^2.*sin(Zeta.*(Y-s)).*exp(-1i*Lambda.*X).*exp(1i*PM*M^2*X);    
else
    
AdRNum = sgn*Zeta.*cos(Zeta.*(Y-s)).*exp(-1i*Lambda.*X);

    if strcmp(type,'pressure')
        AdRNum = -1i*(PM - Lambda).*exp(1i*PM*M^2*X).*AdRNum;
    elseif strcmp(type,'hvelocity')
        AdRNum  = -1i*Lambda.*exp(1i*PM*M^2*X).*AdRNum;
    end
    
end

AdR = AdRNum./AdRDen;


end