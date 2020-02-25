function phi = intField(Z,data,type)

s=data.spac(1); d=data.spac(2); del=data.spac(3);
sigma=data.sigma;
X = real(Z); Y = imag(Z);
TMd=data.TMd;
TPd = data.TPd;
PM0=data.PM0;
GM0 = data.GM0;
KPGM0=data.KPGM0;

omega = data.omega; 
w=data.w;
zetaGM0=mysqrt(omega.*w,GM0);
psiZetaGM0=bsxfun(@times,Y,zetaGM0);
psiMinSZetaGM0=bsxfun(@times,Y-s,zetaGM0);

A1bGM0 = bsxfun(@times,exp(1i*(d*GM0+sigma))./(cos(d*GM0+sigma)-cos(s*zetaGM0)),cos(psiZetaGM0).*exp(-1i*X*GM0));
A1aGM0 = bsxfun(@times,1./(-cos(d*GM0+sigma)+cos(s*zetaGM0)),cos(psiMinSZetaGM0).*exp(-1i*X*GM0));

zeta = @(zVar) mysqrt(omega*w,zVar);

A1aTM = -1./(cos(d*TMd+sigma)-cos(zeta(TMd).*s)).*cos(zeta(TMd).*(Y-s)).*exp(-1i*X.*TMd);
A1bTM = exp(1i*(d*TMd+sigma))./(cos(d*TMd+sigma)-cos(zeta(TMd).*s)).*cos(zeta(TMd).*Y).*exp(-1i*X.*TMd);

A1aTP = -1./(cos(d*TPd+sigma)-cos(zeta(TPd).*s)).*cos(zeta(TPd).*(Y-s)).*exp(-1i*X.*TPd);
A1bTP = exp(1i*(d*TPd+sigma))./(cos(d*TPd+sigma)-cos(zeta(TPd).*s)).*cos(zeta(TPd).*Y).*exp(-1i*X.*TPd);

size(A1aTM),size(data.D1.A)
Aterms = -pi*sum(((data.D1.A+data.D1.C)./(TMd-PM0)).*(A1aTM+A1bTM),3);
Bterms = pi*1i*sum((data.D1.B).*(A1aTP+A1bTP),3);
Tterms = -pi*1i*sum((data.D1.T)./KPGM0.*(A1aGM0+A1bGM0),3);

size(Aterms),size(Bterms),size(Tterms)
phi = Aterms + Bterms + Tterms;

end