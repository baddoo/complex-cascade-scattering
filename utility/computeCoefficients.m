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

TMminGM0 = bsxfun(@minus,permute(TMd,[1,2,4,3,5]),GM0);
coefdata.TMminGM0 = TMminGM0;

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

nfreq=length(kx);
invmat=inv(eye(size(F1))-F1*L1); %this matrix is inv*F

w0  = AAData.Amp(2);
T  = w0/(4i*pi^2*KMGM0);

S = 0;

V = w0*(PM0-GM0)./(4*pi^2*KPGM0);

Tsum = permute(sum(T./TMminGM0,3),[1,2,4,3,5]);

A = (1i*(TMd-PM0)./KpprTM).*Tsum;
A1 = permute(A,[3,4,5,1,2]);
V1 = permute(V,[3,4,5,1,2]);

B1=zeros(size(A1));
C1=zeros(size(A1));
Dres = zeros(size(A1));

% Can remove this looping since it is now obselete
for ifreq = 1:nfreq
Dres(:,:,ifreq) = F1(:,:,ifreq)*A1(:,:,ifreq) + Vpreq(:,:,ifreq)*V1;
B1(:,:,ifreq)=invmat(:,:,ifreq)*Dres(:,:,ifreq);
C1(:,:,ifreq)=L1(:,:,ifreq)*B1(:,:,ifreq);
end

B=permute(B1,[5,4,1,2,3]);
C=permute(C1,[5,4,1,2,3]);

data.T =  T;
data.S =  S;
data.A  = A;
data.B  = B;
data.C  = C;


%Ddata=computeD1Coefficients(coefdata);

data.kx=kx;
data.omega=omega;
data.M=ADData.M;
data.AAData=AAData;
data.ADData=ADData;
data.coefdata = coefdata;
end