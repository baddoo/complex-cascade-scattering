function acousticFieldDownstream = acousticFieldDownstream(data,type,meanFlowData)

LM1 = permute(data.LMa,[3,2,1,4,5]);

A1aResM = data.A1aResM.(type);
A1bResM = data.A1bResM.(type);

TMd=data.TMd;
GM0=data.GM0;
PM0=data.PM0;

%GMm=data.GMm;

KMTM=data.KMTM;
% A1aGM0=data.A1aGM0.(type);
% A1bGM0=data.A1bGM0.(type);

A1aPM0=data.A1aPM0.(type);
A1bPM0=data.A1bPM0.(type);

%KPPM0=data.KPPM0;
KMPM0=data.KMPM0;
KPGM0=data.KPGM0;
KMGM0=data.KMGM0;

%KMGM0=KMGM(:,:,ceil(end/2),:,:);
data.comb=[0,1,0,1];

Dfin= pi*1i*bsxfun(@times,permute(D(LM1,data),[3,2,1,4,5]),(A1aResM + A1bResM));
Dterms = sum(Dfin,3);

gModeFun = A1aPM0+A1bPM0;%A1aGM(:,:,ceil(end/2),:,:)+A1bGM(:,:,ceil(end/2),:,:);

AsumG = -sum((data.D1.A+data.D1.C)./(1i*(TMd-PM0)).*KMTM./KMPM0.*exp(-2i*(TMd-PM0)),3);

%ScoefsG = -(data.(sol).S)./(1i*(GMm-GM0)).*KMGM./KMGM0.*exp(-1i*(GMm-GM0));
%Tterm = -data.D1.T./KPPM0;
%ScoefsG(:,:,ceil(end/2),:,:)=Tterm(:,:,ceil(end/2),:,:);

%SsumG = sum(ScoefsG,3);
TtermsG = -data.D1.T.*KMGM0./KPGM0./KMPM0.*exp(-2i*(GM0-PM0));

P =  AsumG + TtermsG ;

gModeTerms = 1i*pi*P.*gModeFun;
%% Final combination
acousticFieldDownstream = Dterms + gModeTerms ;

end