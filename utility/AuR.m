function AuR = AuR(X,Y,Lambda,Zeta,SQRT,spac,PM,M,type)

s = spac(1); d = spac(2); % Spacing
del = sqrt(s^2 + d^2);
sgn = sign(imag(Lambda(end)));

AuRDen = 2*del*sin(s*Zeta).*SQRT;

if strcmp(type,'vvelocity')
AuRNum =  sgn*Zeta.^2.*exp(1i*sgn*s*Zeta).*sin(Zeta.*Y).*exp(-1i*Lambda.*X).*exp(1i*PM*M^2*X);    
else
    
AuRNum = -sgn*Zeta.*exp(1i*sgn*s*Zeta).*cos(Zeta.*Y).*exp(-1i*Lambda.*X);

    if strcmp(type,'pressure')
        AuRNum = 1i*(PM - Lambda).*exp(1i*PM*M^2*X).*AuRNum;
    elseif strcmp(type,'hvelocity')
        AuRNum  = -1i*Lambda.*exp(1i*PM*M^2*X).*AuRNum;
    end
    
end

AuR = AuRNum./AuRDen;


end