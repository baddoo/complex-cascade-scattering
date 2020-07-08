function [data] = computeCoefficients(ADData,AAData,out)

%% Extract data.
data=out;
KMTM=out.KMTM; KPTP=out.KPTP;
TMd=out.TMd; TPd = out.TPd;
KMGM0= out.KMGM0;
KPGM0= out.KPGM0;
GM0=out.GM0; 
PM0=out.PM0; 
w0  = AAData.Amp(2);
KpprTM=out.KpprTM; KmprTP=out.KmprTP;

L = 1i*(TMd(:) - PM0)./(TPd(:).' - TMd(:)).*KPTP(:).'./KpprTM(:);
F = 1i*exp(2i*(TPd(:) - TMd(:).'))./(TPd(:) - PM0)./(TPd(:) - TMd(:).').*KMTM(:).'./KmprTP(:);

A = (TMd-PM0)*w0./(4*pi^2*KMGM0*KpprTM.*(TMd - GM0));

Dres = F*A(:) ...
       + w0*(PM0-GM0)*exp(2i*(TPd(:) - GM0))./(4i*pi^2*KPGM0*(TPd(:) - GM0).*(TPd(:) - PM0).*KmprTP(:));
B1 = (eye(size(F))-F*L)\Dres;
C1=L*B1;

B=permute(B1,[2,3,1]);
C=permute(C1,[2,3,1]);

data.A  = A;
data.B  = B;
data.C  = C;

data.AAData=AAData;
data.ADData=ADData;
data.comb = out.comb;
end