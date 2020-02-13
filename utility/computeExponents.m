function [newdata,X,Y] = computeExponents(data,plotdata)

newdata=data;

s=data.spac(1); d=data.spac(2); del=data.spac(3);
sigma=data.sigma;
SQRTa=data.SQRTa;
M=data.M;

LMa = data.LMa;
LPa = data.LPa;
ZMa = data.ZMa;
ZPa = data.ZPa;
TMd = data.TMd;
TPd = data.TPd;
GM0 = data.GM0;
PM0 = data.PM0;

omega5 = data.omega5; w=data.w;

X = plotdata.X;
Y = plotdata.Y;
newdata.X = X;
newdata.Y = Y;
% Bit for source modes

xLM3 = bsxfun(@times,X,LMa);         newdata.xLM3 = xLM3;
xLP3 = bsxfun(@times,X,LPa);         newdata.xLP3 = xLP3;
yZM3 = bsxfun(@times,Y,ZMa);         newdata.yZM3 = yZM3;
yZP3 = bsxfun(@times,Y,ZPa);         newdata.yZP3 = yZP3;
yMinhZP3 = bsxfun(@times,Y-s,ZPa);  newdata.yMinhZP3 = yMinhZP3;
yMinhZM3 = bsxfun(@times,Y-s,ZMa);  newdata.yMinhZM3 = yMinhZM3;   

xTM3  = bsxfun(@times,X,TMd);
xTP3  = bsxfun(@times,X,TPd);

phiGM0 = bsxfun(@times,X,GM0);
phiPM0 = bsxfun(@times,X,PM0);

pressureFactor = -1i*exp(1i*X.*PM0.*M^2);
velocityFactor =     exp(1i*X.*PM0.*M^2);

nd3=permute(0:(size(TMd,3)-1),[1,3,2,4,5]);

A1aResP = -bsxfun(@times,ZPa.*exp(1i*s*ZPa)./(del*SQRTa.*sin(s*ZPa)),cos(yZP3).*exp(-1i*xLP3));
newdata.A1aResP.acoustic = A1aResP;
newdata.A1aResP.pressure = A1aResP.*(LPa-PM0).*pressureFactor;
newdata.A1aResP.hvelocity = A1aResP.*(-1i*LPa).*velocityFactor;
newdata.A1aResP.vvelocity = -bsxfun(@times,ZPa.*exp(1i*s*ZPa)./(del*SQRTa.*sin(s*ZPa)).*ZPa,-sin(yZP3).*exp(-1i*xLP3)).*velocityFactor;

A1bResP = bsxfun(@times,ZPa./(del*SQRTa.*sin(s*ZPa)),cos(yMinhZP3).*exp(-1i*xLP3));
newdata.A1bResP.acoustic = A1bResP;
newdata.A1bResP.pressure = A1bResP.*(LPa-PM0).*pressureFactor;
newdata.A1bResP.hvelocity = A1bResP.*(-1i*LPa).*velocityFactor;
newdata.A1bResP.vvelocity = bsxfun(@times,ZPa./(del*SQRTa.*sin(s*ZPa)).*ZPa,-sin(yMinhZP3).*exp(-1i*xLP3)).*velocityFactor;

A1aResM = -bsxfun(@times,ZMa.*exp(-1i*s*ZMa)./(del*SQRTa.*sin(s*ZMa)),cos(yZM3).*exp(-1i*xLM3));
newdata.A1aResM.acoustic=A1aResM;
newdata.A1aResM.pressure=A1aResM.*(LMa-PM0).*pressureFactor;
newdata.A1aResM.hvelocity=A1aResM.*(-1i*LMa).*velocityFactor;
newdata.A1aResM.vvelocity=-bsxfun(@times,ZMa.*exp(-1i*s*ZMa)./(del*SQRTa.*sin(s*ZMa)).*ZMa,-sin(yZM3).*exp(-1i*xLM3)).*velocityFactor;

A1bResM = bsxfun(@times,ZMa./(del*SQRTa.*sin(s*ZMa)),cos(yMinhZM3).*exp(-1i*xLM3));
newdata.A1bResM.acoustic = A1bResM;
newdata.A1bResM.pressure = A1bResM.*(LMa-PM0).*pressureFactor;
newdata.A1bResM.hvelocity = A1bResM.*(-1i*LMa).*velocityFactor;
newdata.A1bResM.vvelocity = bsxfun(@times,ZMa./(del*SQRTa.*sin(s*ZMa)).*ZMa,-sin(yMinhZM3).*exp(-1i*xLM3)).*velocityFactor;

zeta = @(zVar) mysqrt(omega5*w,zVar);

A1bTP = bsxfun(@times,exp(1i*(d*TPd+sigma))./(cos(d*TPd+sigma)-cos(zeta(TPd).*s)),cos(zeta(TPd).*Y).*exp(-1i*xTP3));
newdata.A1bTP.acoustic = A1bTP;
newdata.A1bTP.pressure = A1bTP.*(TPd-PM0).*pressureFactor;
newdata.A1bTP.hvelocity = A1bTP.*(-1i*TPd).*velocityFactor;
newdata.A1bTP.vvelocity =  bsxfun(@times,exp(1i*(d*TPd+sigma))./(cos(d*TPd+sigma)-cos(zeta(TPd).*s)),-zeta(TPd).*sin(zeta(TPd).*Y).*exp(-1i*xTP3)).*velocityFactor;

A1aTP = -bsxfun(@times,1./(cos(d*TPd+sigma)-cos(zeta(TPd).*s)),cos(zeta(TPd).*(Y-s)).*exp(-1i*xTP3));
newdata.A1aTP.acoustic = A1aTP;
newdata.A1aTP.pressure = A1aTP.*(TPd-PM0).*pressureFactor;
newdata.A1aTP.hvelocity = A1aTP.*(-1i*TPd).*velocityFactor;
newdata.A1aTP.vvelocity =  -bsxfun(@times,1./(cos(d*TPd+sigma)-cos(zeta(TPd).*s)),-zeta(TPd).*sin(zeta(TPd).*(Y-s)).*exp(-1i*xTP3)).*velocityFactor;

A1bTM = bsxfun(@times,exp(1i*(d*TMd+sigma))./(cos(d*TMd+sigma)-cos(zeta(TMd).*s)),cos(zeta(TMd).*Y).*exp(-1i*xTM3));
newdata.A1bTM.acoustic = A1bTM;
newdata.A1bTM.pressure = A1bTM.*(TMd-PM0).*pressureFactor;
newdata.A1bTM.hvelocity = A1bTM.*(-1i*TMd).*velocityFactor;
newdata.A1bTM.vvelocity = bsxfun(@times,exp(1i*(d*TMd+sigma))./(cos(d*TMd+sigma)-cos(zeta(TMd).*s)),-zeta(TMd).*sin(zeta(TMd).*Y).*exp(-1i*xTM3)).*velocityFactor;

A1aTM = -bsxfun(@times,1./(cos(d*TMd+sigma)-cos(zeta(TMd).*s)),cos(zeta(TMd).*(Y-s)).*exp(-1i*xTM3));
newdata.A1aTM.acoustic = A1aTM;
newdata.A1aTM.pressure = A1aTM.*(TMd-PM0).*pressureFactor;
newdata.A1aTM.hvelocity = A1aTM.*(-1i*TMd).*velocityFactor;
newdata.A1aTM.vvelocity = -bsxfun(@times,1./(cos(d*TMd+sigma)-cos(zeta(TMd).*s)),-zeta(TMd).*sin(zeta(TMd).*(Y-s)).*exp(-1i*xTM3)).*velocityFactor;

zetaGM0=mysqrt(omega5.*w,GM0);
psiZetaGM0=bsxfun(@times,Y,zetaGM0);
psiMinhZetaGM0=bsxfun(@times,Y-s,zetaGM0);

A1bGM0 = bsxfun(@times,exp(1i*(d*GM0+sigma))./(cos(d*GM0+sigma)-cos(s*zetaGM0)),cos(psiZetaGM0).*exp(-1i*phiGM0));
newdata.A1bGM0.acoustic = A1bGM0;
newdata.A1bGM0.pressure = A1bGM0.*(GM0-PM0).*pressureFactor;
newdata.A1bGM0.hvelocity = A1bGM0.*(-1i*GM0).*velocityFactor;
newdata.A1bGM0.vvelocity = bsxfun(@times,exp(1i*(d*GM0+sigma))./(cos(d*GM0+sigma)-cos(s*zetaGM0)).*zetaGM0,-sin(psiZetaGM0).*exp(-1i*phiGM0)).*velocityFactor;
 
A1aGM0 = bsxfun(@times,1./(-cos(d*GM0+sigma)+cos(s*zetaGM0)),cos(psiMinhZetaGM0).*exp(-1i*phiGM0));
newdata.A1aGM0.acoustic = A1aGM0;
newdata.A1aGM0.pressure = A1aGM0.*(GM0-PM0).*pressureFactor;
newdata.A1aGM0.hvelocity = A1aGM0.*(-1i*GM0).*velocityFactor;
newdata.A1aGM0.vvelocity = bsxfun(@times,-1./(cos(d*GM0+sigma)-cos(s*zetaGM0)).*zetaGM0,-sin(psiMinhZetaGM0).*exp(-1i*phiGM0)).*velocityFactor;

zetaPM0=mysqrt(omega5.*w,PM0);
psiZetaPM0=bsxfun(@times,Y,zetaPM0);
psiMinhZetaPM0=bsxfun(@times,Y-s,zetaPM0);
 
A1bPM0 = bsxfun(@times,exp(1i*(d*PM0+sigma))./(cos(d*PM0+sigma)-cos(s*zetaPM0)),cos(psiZetaPM0).*exp(-1i*phiPM0));
newdata.A1bPM0.acoustic = A1bPM0;
newdata.A1bPM0.pressure = 0;
newdata.A1bPM0.hvelocity = A1bPM0.*(-1i*PM0).*velocityFactor;
newdata.A1bPM0.vvelocity = bsxfun(@times,exp(1i*(d*PM0+sigma))./(cos(d*PM0+sigma)-cos(s*zetaPM0)).*zetaPM0,-sin(psiZetaPM0).*exp(-1i*phiPM0)).*velocityFactor;
 
A1aPM0 = bsxfun(@times,1./(-cos(d*PM0+sigma)+cos(s*zetaPM0)),cos(psiMinhZetaPM0).*exp(-1i*phiPM0));
newdata.A1aPM0.acoustic = A1aPM0;
newdata.A1aPM0.pressure = 0;
newdata.A1aPM0.hvelocity = A1aPM0.*(-1i*PM0).*velocityFactor;
newdata.A1aPM0.vvelocity = bsxfun(@times,-1./(cos(d*PM0+sigma)-cos(s*zetaPM0)).*zetaPM0,-sin(psiMinhZetaPM0).*exp(-1i*phiPM0)).*velocityFactor;

end
