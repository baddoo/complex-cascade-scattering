function [KpprTM,KmprTP]=computeKDerivatives(ADData,AAData,TMd,KMTM,TPd,KPTP)
  
s=ADData.spac(1); d=ADData.spac(2);
Sigma = AAData.Sigma;
omega = AAData.omega;

w = AAData.w;
mu = ADData.mu;

zTM = mysqrt(omega*w,TMd);
zTP = mysqrt(omega*w,TPd);

TMden = cos(s*zTM)-cos(d*TMd+Sigma);
TMnum = (-TMd./zTM.*sin(s*zTM)-s*TMd.*cos(s*zTM)).*TMden...
        -zTM.*sin(s.*zTM).*(TMd.*s./zTM.*sin(s*zTM)+d*sin(d*TMd+Sigma));
KprTM=  (TMnum./TMden.^2 - 1i*mu(2) - 2*mu(3)*TMd)/4/pi;

KpprTM=KprTM./KMTM;

TPden = cos(s*zTP)-cos(d*TPd+Sigma);
TPnum = (-TPd./zTP.*sin(s*zTP)-s*TPd.*cos(s*zTP)).*TPden...
        -zTP.*sin(s.*zTP).*(TPd.*s./zTP.*sin(s*zTP)+d*sin(d*TPd+Sigma));
KprTP=  (TPnum./TPden.^2 - 1i*mu(2) - 2*mu(3)*TPd )/4/pi;

KmprTP=KprTP./KPTP;

end