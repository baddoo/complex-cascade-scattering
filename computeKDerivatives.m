function [KpprTM,KmprTP]=computeKDerivatives(ADData,AAData,TMd,KMTM,TPd,KPTP)
  
he=ADData.spac(1); d=ADData.spac(2);
sigma5 = AAData.sigma5;
omega5 = AAData.omega5;

w = AAData.w;
mu = ADData.mu;

zTM = mysqrt(omega5*w,TMd);
zTP = mysqrt(omega5*w,TPd);

TMden = cos(he*zTM)-cos(d*TMd+sigma5);
TMnum = (-TMd./zTM.*sin(he*zTM)-he*TMd.*cos(he*zTM)).*TMden...
        -zTM.*sin(he.*zTM).*(TMd.*he./zTM.*sin(he*zTM)+d*sin(d*TMd+sigma5));
KprTM=  (TMnum./TMden.^2 -1i*mu(2) - 2*mu(3)*TMd)/4/pi;

KpprTM=KprTM./KMTM;

TPden = cos(he*zTP)-cos(d*TPd+sigma5);
TPnum = (-TPd./zTP.*sin(he*zTP)-he*TPd.*cos(he*zTP)).*TPden...
        -zTP.*sin(he.*zTP).*(TPd.*he./zTP.*sin(he*zTP)+d*sin(d*TPd+sigma5));
KprTP=  (TPnum./TPden.^2 -1i*mu(2) - 2*mu(3)*TPd )/4/pi;

KmprTP=KprTP./KPTP;

end