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