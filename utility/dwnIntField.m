function phi = dwnIntField(Z,data,type)

s=data.spac(1); d=data.spac(2); del=data.spac(3);
sigma=data.sigma;
omega = data.omega; 
w=data.w;
SQRTa=data.SQRTa;

X = real(Z); Y = imag(Z);

TMd=data.TMd;
TPd=data.TPd;

GM0=data.GM0;
PM0=data.PM0;

KMTM=data.KMTM;
% 
KPGM0=data.KPGM0;
KMGM0=data.KMGM0;

KMPM0=data.KMPM0;

LMa = data.LMa;
LMa = data.LMa;

ZMa = data.ZMa;


data.comb=[0,1,0,1];
 
yMinSZM3 = bsxfun(@times,Y-s,ZMa); 

Dcoefs = permute(D(permute(LMa,[3,2,1]),data),[3,2,1]);

Dfin = pi*1i*Dcoefs.*ZMa./(del*SQRTa.*sin(s*ZMa)).*cos(yMinSZM3).*exp(-1i*X.*LMa);
Dterms = sum(Dfin,3);

A1bTM = exp(1i*(d*TMd+sigma))./(cos(d*TMd+sigma)-cos(zeta(TMd).*s)).*cos(zeta(TMd).*Y).*exp(-1i*X.*TMd);
Aterms = -pi*sum((data.D1.A+data.D1.C)./(TMd-PM0).*A1bTM,3);

A1bTP = bsxfun(@times,exp(1i*(d*TPd+sigma))./(cos(d*TPd+sigma)-cos(zeta(TPd).*s)),cos(zeta(TPd).*Y).*exp(-1i*X.*TPd));
Bterms = pi*1i*sum((data.D1.B).*A1bTP,3);
 
AsumG = -sum((data.D1.A+data.D1.C)./(1i*(TMd-PM0)).*KMTM./KMPM0.*exp(-2i*(TMd-PM0)),3);

zetaGM0=mysqrt(omega.*w,GM0);
psiZetaGM0=bsxfun(@times,Y,zetaGM0);
A1bGM0 = bsxfun(@times,exp(1i*(d*GM0+sigma))./(cos(d*GM0+sigma)-cos(s*zetaGM0)),cos(psiZetaGM0).*exp(-1i*X*GM0));
Tsum2 = data.D1.T./KPGM0.*A1bGM0; 
Tterms2 = -pi*1i*sum(Tsum2,3);

TtermsG = -data.D1.T.*KMGM0./KPGM0./KMPM0.*exp(-2i*(GM0-PM0));

P =  AsumG + TtermsG ;

zetaPM0=mysqrt(omega.*w,PM0);
yMinSZetaPM0=bsxfun(@times,Y-s,zetaPM0);
A1aPM0 = bsxfun(@times,1./(-cos(d*PM0+sigma)+cos(s*zetaPM0)),cos(yMinSZetaPM0).*exp(-1i*X*PM0));

gModeTerms = 1i*pi*P.*A1aPM0;
%% Final sums
phi = Dterms;% + 0*Aterms + 0*Bterms + 0*Tterms2  + 0*gModeTerms;

end