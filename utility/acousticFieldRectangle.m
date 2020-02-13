function acousticFieldRectangle = acousticFieldRectangle(data,type)

A1aTM=data.A1aTM.(type);
A1aTP=data.A1aTP.(type);
A1bTM=data.A1bTM.(type);
A1bTP=data.A1bTP.(type);
TMd=data.TMd;
GM0=data.GM0;
PM0=data.PM0;
A1aGM0=data.A1aGM0.(type);
A1bGM0=data.A1bGM0.(type);
KPGM0=data.KPGM0;
%KPGM0=data.KPGM0;

Aterms = -pi*sum(((data.D1.A+data.D1.C)./(TMd-PM0)).*(A1aTM+A1bTM),3);
Bterms = pi*1i*sum((data.D1.B).*(A1aTP+A1bTP),3);
Tterms = -pi*1i*sum((data.D1.T)./KPGM0.*(A1aGM0+A1bGM0),3);
%Vterms = 0*pi*1i*data.D1.V.*data.KPPM0.*data.A1bPM0.(type);

acousticFieldRectangle= Aterms + Bterms + Tterms;

end