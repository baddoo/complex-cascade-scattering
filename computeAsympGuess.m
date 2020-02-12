function asympGuess = computeAsympGuess(ADData,AAData,Modes)

mu = ADData.mu;
he = ADData.spac(1);
d = ADData.spac(2);
s = sqrt(he^2+d^2);
sigma5 = AAData.sigma5;

n= 1:Modes.trunc;

if mu(2) ==0 && mu(3) == 0
% No-mean-flow Darcy
    a  = (d+1i*he)*2*pi/s.^2;
    b  = 1i*(d+1i*he)/s.^2;
    c  = (pi-sigma5+1i*log(2*pi/mu(1)/s)-atan(he/d))*(d+1i*he)./s.^2;
    c4 = conj(c);
    
    a1 = (-d+1i*he)*2*pi/s.^2;
    b1 =  -1i*(-d+1i*he)/s.^2;
    c1 = (-pi+sigma5-1i*log(2*pi/mu(1)/s)-atan(he/d))*(-d+1i*he)./s.^2;
    c3 = conj(c1);
    
elseif mu(3)==0
% Mean-flow Darcy
    a  = (d+1i*he)*2*pi/s.^2;
    b  = 0;
    c = ( -sigma5+1i*log(1+1./(1i*mu(2))))*(d+1i*he)./s.^2;

    a1 = (-d+1i*he)*2*pi/s.^2;
    b1 =  0;
    c1 = ( sigma5-1i*log(1-1./(1i*mu(2))))*(-d+1i*he)./s.^2;
    
    c4 = (-sigma5-1i*log(1+1./(1i*mu(2))))*(d-1i*he)./s.^2;
    c3 = ( sigma5+1i*log(1-1./(1i*mu(2))))*(-d-1i*he)./s.^2;

else
% Impedance 
    a  = (d+1i*he)*2*pi/s.^2;
    b  = 0;
    c  = -sigma5*(d+1i*he)./s.^2;

    a1 = (-d+1i*he)*2*pi/s.^2;
    b1 =  0;
    c1 = sigma5*(-d+1i*he)./s.^2;
    
    c4 = -sigma5*( d-1i*he)./s.^2;
    c3 =  sigma5*(-d-1i*he)./s.^2;

end

Rsols1b = a*n+b*log(n)+c;
Rsols2b = a1*n+b1*log(n)+c1;
Rsols3b = conj(a1)*n+conj(b1)*log(n)+c3;
Rsols4b = conj(a)*n+conj(b)*log(n)+c4;
asympGuess = [Rsols1b,Rsols2b,Rsols3b,Rsols4b];

end