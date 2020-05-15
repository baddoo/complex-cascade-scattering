function D = D(gam,data)

newdata=data;

newdata.gam=gam;

Kargs=data.Kargs;
newdata.gminGM0 = bsxfun(@minus,gam,data.GM0);
newdata.gminPM0 = bsxfun(@minus,gam,data.PM0);

newdata.gminTM = bsxfun(@minus,gam,data.TMd);
newdata.gminTP = bsxfun(@minus,gam,data.TPd);

newdata.gminLM = bsxfun(@minus,gam,data.LMa);
newdata.gminLP = bsxfun(@minus,gam,data.LPa);

newdata.KPG = Kplus(gam,Kargs);
newdata.KMG = Kminus(gam,Kargs);


D=combineD(newdata,newdata.D1);

function D1 = D1(data,data2)

gminGM0 = data.gminGM0;

KPG=data.KPG;
T = data2.T;

Tsum=sum(T./gminGM0,3)./KPG;

D1 = Tsum;

end

function D3=D3(data,data2)

B=data2.B;
gminTP=data.gminTP;
KPG=data.KPG;
KPTP=data.KPTP;

num=B.*KPTP;
denom=gminTP;

D3terms=bsxfun(@rdivide,num,denom);

D3=-sum(D3terms,3)./KPG;

end

function [D2] = D2(data,data2)

A=data2.A;
V=data2.V;

%S=data2.S;

KMTM=data.KMTM;
%KMGM0=data.KMGM0;
KMG=data.KMG;
gminGM0 = data.gminGM0;
gminPM0 = data.gminPM0;
gminTM = data.gminTM;

D2V = V.*exp(2i*gminGM0)./(1i*gminGM0.*gminPM0.*KMG);

D2termsA=bsxfun(@rdivide,A.*KMTM,gminTM).*exp(2i*gminTM);
D2A=-sum(D2termsA,3)./(1i*KMG.*gminPM0);

D2 = D2V + D2A;

end

function [D4]=D4(data,data2)

C=data2.C;
KMG=data.KMG;
gminPM0=data.gminPM0;
KMTM=data.KMTM;

gminTM=data.gminTM;

D4terms=bsxfun(@rdivide,C.*KMTM,gminTM).*exp(2i*gminTM);
D4=-sum(D4terms,3)./(1i*KMG.*gminPM0);

end

function D=combineD(data,dnum)

if data.comb(1)==1; D_O1 = D1(data,dnum); else D_O1=0; end %#ok<*SEPEX>
if data.comb(2)==1; D_O2 = D2(data,dnum); else D_O2=0; end
if data.comb(3)==1; D_O3 = D3(data,dnum); else D_O3=0; end
if data.comb(4)==1; D_O4 = D4(data,dnum); else D_O4=0; end

D=D_O1+D_O2+D_O3+D_O4;

end

end

