function [data] = computeCoefficients(ADData,AAData,out)

%% Extract data.
data=out;
k5=AAData.k5; omega5=AAData.omega5;
KMTM=out.KMTM; KPTP=out.KPTP;
TMd=out.TMd; TPd = out.TPd;
LMa=out.LMa; LPa = out.LPa; 

dmodes=out.dmodes;
GM0=out.GM0; 
PM0=out.PM0; 

KpprTM=out.KpprTM; KmprTP=out.KmprTP;

coefdata=out;
coefdata.AAData = AAData;
%coefdata.Amp5=AAData.Amp5;
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
F=-exp(1i*permute(TPminTM,[1,2,4,3,5]))./(1i*permute(bsxfun(@times,TPd-PM0,TPminTM),[1,2,4,3,5])).*finmatF;

finmatV=bsxfun(@rdivide,1,permute(KmprTP,[1,2,4,3,5]));
Vpre=exp(1i*TPminGM0)./(1i*permute(bsxfun(@times,(TPd-GM0),(TPd-PM0)),[1,2,4,3,5])).*finmatV;

Vpreq= permute(Vpre,[4,3,5,1,2]);
F1= permute(F,[4,3,5,1,2]);
L1= permute(L,[4,3,5,1,2]);

coefdata.F1= F1;
coefdata.Vpreq= Vpreq;
coefdata.L1= L1;

nfreq=length(k5);
invmat=zeros(size(F1));

for ifreq = 1:nfreq
invmat(:,:,ifreq)=inv(eye(dmodes)-F1(:,:,ifreq)*L1(:,:,ifreq)); %this matrix is inv*F
end
coefdata.invmat=invmat;

data.D1=computeD1Coefficients(coefdata);

data.k5=k5;
data.omega5=omega5;
data.delta=ADData.delta;
data.M=ADData.M;
%data.Amp5=AAData.Amp5;
data.AAData=AAData;
data.ADData=ADData;
data.coefdata = coefdata;
end