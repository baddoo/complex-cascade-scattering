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

end