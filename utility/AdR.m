function AdR = AdR(X,Y,Lambda,Zeta,SQRT,spac)

s = spac(1); d = spac(2); % Spacing
del = sqrt(s^2 + d^2);
sgn = sign(imag(Lambda(end)));

AdRNum = sgn*Zeta.*cos(Zeta.*(Y-s)).*exp(-1i*Lambda.*X);
AdRDen = 2*del*sin(s*Zeta).*SQRT;

AdR = AdRNum./AdRDen;

end