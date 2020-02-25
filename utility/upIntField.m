function phi = upIntField(Z,data,type)

X = real(Z); Y = imag(Z);
s=data.spac(1); d=data.spac(2); del=data.spac(3);
omega = data.omega; 
w=data.w;

sigma=data.sigma;

LPa = data.LPa;
ZPa = data.ZPa;
TMd = data.TMd;
TPd = data.TPd;

SQRTa=data.SQRTa;

GM0 = data.GM0;
zetaGM0=mysqrt(omega.*w,GM0);
zeta = @(zVar) mysqrt(omega*w,zVar);

PM0=data.PM0;
A1aTM=-bsxfun(@times,1./(cos(d*TMd+sigma)-cos(zeta(TMd).*s)),cos(zeta(TMd).*(Y-s)).*exp(-1i*X.*TMd));
A1aTP=-bsxfun(@times,1./(cos(d*TPd+sigma)-cos(zeta(TPd).*s)),cos(zeta(TPd).*(Y-s)).*exp(-1i*X.*TPd));

%A1aGM=data.A1aGM.(type);
%KPGM=data.KPGM;

data.comb=[1,0,1,0];
%A1aResP = data.A1aResP.(type);
%Dfin= pi*1i*bsxfun(@times,permute(D(LPa,data),[3,2,1,4,5]),(A1aResP));

AplusC = (data.D1.A+data.D1.C);
Aterms = -pi*sum(AplusC./(TMd-PM0).*A1aTM,3);
%Aterms = -pi*sum((data.D1.A)./(TMd-GM0).*(A1aTM+0*A1bTM),3);

B = data.D1.B;
Bterms = pi*1i*sum(B.*A1aTP,3);
T = data.D1.T;

psiMinSZetaGM0=bsxfun(@times,Y-s,zetaGM0);

A1aGM0 = bsxfun(@times,1./(-cos(d*GM0+sigma)+cos(s*zetaGM0)),cos(psiMinSZetaGM0).*exp(-1i*X*GM0));
Tterms = -pi*1i*T./data.KPGM0.*A1aGM0;
%Tterms = 0;

%upIntField = sum(Dfin,3)+Aterms+Bterms+Tterms;%+Vterms;

Dcoefs = permute(D(permute(LPa,[3,2,1]),data),[3,2,1]);
%phi = @(zVar) pi/del*sum(permute(D(LPa,data),[3,2,1])./SQRTa.*exp(-1i*(X*LPa-Y*ZPa)),3);
phi = 2*pi/del.*sum(Dcoefs.*ZPa./SQRTa./(1-exp(-2i*s*ZPa)).*cos(Y.*ZPa).*exp(-1i*X.*LPa),3)...   
      + Aterms + Bterms + Tterms;

end