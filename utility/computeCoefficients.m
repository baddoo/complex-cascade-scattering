function [data] = computeCoefficients(ADData,AAData,out)

%% Extract data.
data=out;
kx=AAData.kx; omega=AAData.omega;
KMTM=out.KMTM; KPTP=out.KPTP;
TMd=out.TMd; TPd = out.TPd;
LMa=out.LMa; LPa = out.LPa; 
KMGM0= out.KMGM0;
KPGM0= out.KPGM0;
GM0=out.GM0; 
PM0=out.PM0; 

KpprTM=out.KpprTM; KmprTP=out.KmprTP;

coefdata=out;
coefdata.AAData = AAData;
coefdata.Beta=ADData.Beta;

%% Define differences between modes

TPminTM=bsxfun(@minus,TPd,permute(TMd,[1,2,4,3,5]));
coefdata.TPminTM=TPminTM;

TPminGM0 = bsxfun(@minus,permute(TPd,[1,2,4,3,5]),GM0);
coefdata.TPminGM0 = TPminGM0;

TMminLM = bsxfun(@minus,permute(TMd,[1,2,4,3,5]),LMa);
coefdata.TMminLM = TMminLM;

coefdata.TPminLP=bsxfun(@minus,permute(TPd,[1,2,4,3,5]),LPa);

finmatL=bsxfun(@rdivide,KPTP,permute(KpprTM,[1,2,4,3,5])); 
L= 1i*bsxfun(@rdivide,permute(TMd-PM0,[1,2,4,3,5]),TPminTM).*finmatL;

finmatF=bsxfun(@rdivide,KMTM,permute(KmprTP,[1,2,4,3,5]));
F=-exp(2i*permute(TPminTM,[1,2,4,3,5]))./(1i*permute(bsxfun(@times,TPd-PM0,TPminTM),[1,2,4,3,5])).*finmatF;

finmatV=bsxfun(@rdivide,1,permute(KmprTP,[1,2,4,3,5]));
Vpre=exp(2i*TPminGM0)./(1i*permute(bsxfun(@times,(TPd-GM0),(TPd-PM0)),[1,2,4,3,5])).*finmatV;

Vpreq= permute(Vpre,[4,3,5,1,2]);
F1= permute(F,[4,3,5,1,2]);
L1= permute(L,[4,3,5,1,2]);

w0  = AAData.Amp(2);

V = w0*(PM0-GM0)./(4*pi^2*KPGM0);

A = (1i*(TMd-PM0)./KpprTM).*w0./(4i*pi^2*KMGM0.*(TMd - GM0));
A1 = permute(A,[3,4,5,1,2]);
V1 = permute(V,[3,4,5,1,2]);

Dres = F1*A1 + Vpreq*V1;
B1 = (eye(size(F1))-F1*L1)\Dres;
C1=L1*B1;

B=permute(B1,[5,4,1,2,3]);
C=permute(C1,[5,4,1,2,3]);

data.A  = A;
data.B  = B;
data.C  = C;

data.kx=kx;
data.omega=omega;
data.M=ADData.M;
data.AAData=AAData;
data.ADData=ADData;
data.coefdata = coefdata;
end