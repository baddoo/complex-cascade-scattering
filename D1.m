function D1 = D1(data,data2)

gminGM0 = data.gminGM0;

KPG=data.KPG;
T = data2.T;

Tsum=sum(T./gminGM0,3)./KPG;

D1 = Tsum;

end