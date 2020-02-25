function acousticFieldDownstream = acousticFieldDownstream(data,type,meanFlowData)

LM1 = permute(data.LMa,[3,2,1,4,5]);

A1aResM = data.A1aResM.(type);
A1bResM = data.A1bResM.(type);

TMd=data.TMd;
GM0=data.GM0;
PM0=data.PM0;

KMTM=data.KMTM;

A1aPM0=data.A1aPM0.(type);
A1bPM0=data.A1bPM0.(type);

KMPM0=data.KMPM0;
KPGM0=data.KPGM0;
KMGM0=data.KMGM0;

data.comb=[0,1,0,1];

Dfin= pi*1i*bsxfun(@times,permute(D(LM1,data),[3,2,1,4,5]),(A1aResM + A1bResM));
Dterms = sum(Dfin,3);

gModeFun = A1aPM0 + A1bPM0;

AsumG = -sum((data.D1.A+data.D1.C)./(1i*(TMd-PM0)).*KMTM./KMPM0.*exp(-2i*(TMd-PM0)),3);

TtermsG = -data.D1.T.*KMGM0./KPGM0./KMPM0.*exp(-2i*(GM0-PM0));

P =  AsumG + TtermsG ;

gModeTerms = 1i*pi*P.*gModeFun;

%% Final combination
acousticFieldDownstream = Dterms + gModeTerms ;

end