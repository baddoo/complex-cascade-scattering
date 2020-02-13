function acousticFieldUpperTriangle = acousticFieldUpperTriangle(data,type)

A1aResP = data.A1aResP.(type);

LPa = permute(data.LPa,[3,2,1,4,5]);
TMd = data.TMd;
GM0=data.GM0;
PM0=data.PM0;
A1aTM=data.A1aTM.(type);
A1aTP=data.A1aTP.(type);

%A1aGM=data.A1aGM.(type);
%KPGM=data.KPGM;

data.comb=[1,0,1,0];
Dfin= pi*1i*bsxfun(@times,permute(D(LPa,data),[3,2,1,4,5]),(A1aResP));

Aterms = -pi*sum((data.D1.A+data.D1.C)./(TMd-PM0).*A1aTM,3);
%Aterms = -pi*sum((data.D1.A)./(TMd-GM0).*(A1aTM+0*A1bTM),3);

Bterms = pi*1i*sum(data.D1.B.*A1aTP,3);
Tterms = -pi*1i*data.D1.T./data.KPGM0.*data.A1aGM0.(type);
%Tterms = 0;

acousticFieldUpperTriangle=sum(Dfin,3)+Aterms+Bterms+Tterms;%+Vterms;

end