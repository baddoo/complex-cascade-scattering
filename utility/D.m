function D = D(gam,data)

Beta = data.ADData.Beta;
w0 = data.AAData.Amp(2);
PM0 = data.PM0;
GM0 = data.GM0;

Kargs=data.Kargs;
gminGM0 = gam - GM0;
gminPM0 = gam - PM0;

gminTM = gam - data.TMd;
gminTP = gam - data.TPd;

KPG = Kplus(gam,Kargs);
KMG = Kminus(gam,Kargs);

KPTP = data.KPTP;    
KMTM = data.KMTM;

KPGM0 = data.KPGM0;
KMGM0 = data.KMGM0;

if data.comb(1)==1
    
    w0 = data.AAData.Amp(2);
    D1 = w0./(4i*pi^2*gminGM0.*KMGM0.*KPG);

else
    D1=0;
end

if data.comb(2)==1

A = data.A;

D2termsA = bsxfun(@rdivide,A.*KMTM,gminTM).*exp(2i*gminTM);
D2A = -sum(D2termsA,3)./(1i*KMG.*gminPM0);

D2 = w0*Beta^2*(PM0-GM0)*exp(2i*gminGM0)./(4i*pi^2*Beta^2*KPGM0*gminGM0.*gminPM0.*KMG) ...
    +  D2A;  
    
else
    D2=0;
end

if data.comb(3)==1
    
B=data.B;
D3terms=B.*KPTP./gminTP;
D3 = -sum(D3terms,3)./KPG;

else
    D3=0; 
end

if data.comb(4)==1
    
C=data.C;
D4terms=C.*KMTM./gminTM.*exp(2i*gminTM);
D4 = -sum(D4terms,3)./(1i*KMG.*gminPM0);

else
    D4=0; 
end

D=D1+D2+D3+D4;

end

