function output = computeModes(ADData,AAData,Modes)
% computeModes Computes the relevant modes for the problem

%% Get data from structures
s=ADData.spac(1); d=ADData.spac(2); del = ADData.spac(3);
trunc=Modes.trunc; w=AAData.w; Sigma=AAData.Sigma;
kx=AAData.kx; omega = AAData.omega; chie = ADData.chie;
Beta = ADData.Beta; amodes=Modes.amodes;

%% Compute duct mods
R0 = 5;
err = inf;
counter = 0;

while err>1e-2
chi = 0;
if cos(chi)<0; error('The ellipse does not have the right parameters'); end
[TP,TM,asymP] = findDuctModes(R0,chi,ADData,AAData,Modes);

% % Remove the largest roots so we have as many duct modes as acoustic modes
% if numel(TP)>2*trunc+1
% TP = mink(TP,2*trunc+1,'ComparisonMethod','abs');
% else
% trunc = (min(numel(TP),numel(TM))-1)/2;
% end
% if numel(TM)>2*trunc+1
% TM = mink(TM,2*trunc+1,'ComparisonMethod','abs');
% else
% trunc = (min(numel(TP),numel(TM))-1)/2;
% end
Modes.trunc = trunc;

%% Compute acoustic modes
aTrunc = permute(-trunc:trunc,[1,3,2]);
f=(bsxfun(@minus,Sigma,2*pi*aTrunc))/del;
SQRT=mysqrt(w*omega,f); output.SQRT=SQRT;
LM=-f*sin(chie)-cos(chie)*SQRT;
LP=-f*sin(chie)+cos(chie)*SQRT;

%% Compute other modes
dmodes = Modes.dmodes;

TMd = TM(:,:,(1:dmodes),:,:);
TPd = TP(:,:,(1:dmodes),:,:);

ZM= -(Sigma+d*LM-2*pi*permute(-trunc:trunc,[1,3,2,4,5]))/s;
ZP=  (Sigma+d*LP-2*pi*permute(-trunc:trunc,[1,3,2,4,5]))/s;

LMa= LM(:,:,(ceil(end/2)-amodes):(ceil(end/2)+amodes),:,:); 
LPa= LP(:,:,(ceil(end/2)-amodes):(ceil(end/2)+amodes),:,:); 
ZMa= ZM(:,:,(ceil(end/2)-amodes):(ceil(end/2)+amodes),:,:); 
ZPa= ZP(:,:,(ceil(end/2)-amodes):(ceil(end/2)+amodes),:,:); 

SQRTa = SQRT(:,:,(ceil(end/2)-amodes):(ceil(end/2)+amodes),:,:);

Kargs.ADData=ADData; Kargs.AAData=AAData;
Kargs.TM=TM; Kargs.LM=LM; Kargs.TP=TP; Kargs.LP=LP;

% Check factorisation. If the error is too small, then go back and increase
% the radius of integration by 50%.

largeVal = sqrt(R0)*exp(1i*linspace(-pi,pi,20));
err = norm(Kminus(largeVal,Kargs).*Kplus(largeVal,Kargs) - regK(largeVal,ADData,AAData),'inf');

counter = counter + 1;

R0 = 1.5*R0;

if counter>10; error("The roots couldn't all be found"); end

end

KMTM=permute(Kminus(permute(TMd,[3,2,1,4,5]),Kargs),[3,2,1,4,5]);
KPTP=permute(Kplus(permute(TPd,[3,2,1,4,5]),Kargs),[3,2,1,4,5]); 

KMLM=permute(Kminus(permute(LMa,[3,2,1,4,5]),Kargs),[3,2,1,4,5]);
KPLP=permute(Kplus(permute(LPa,[3,2,1,4,5]),Kargs),[3,2,1,4,5]); 

GM0 = -(kx)/Beta^2;
PM0 = -(omega)/Beta^2;

output.GM0 = GM0;
output.PM0 = PM0;

%% Phase plot for debugging.
pl=0;
if pl == 1
    figure(1)

% Define grid for phase plot
nx = 100; ny = 100;
xp = linspace(-25,15,nx);
yp = -linspace(-15,15,ny);
[X,Y]=meshgrid(xp,yp);
Z = X + 1i*Y;
% Plot phase plot of error. Should be all red, which means that the error
% is zero.
Kerr = regK(Z,ADData,AAData)-(Kplus(Z,Kargs).*Kminus(Z,Kargs));
%PhasePlot(Z,Kerr,'d');
pcolor(real(Z),imag(Z),abs(Kerr));
colormap gray
shading interp

hold on
plot(TM(:),'bx') % Roots in LHP in blue
plot(TP(:),'rx') % Roots in UHP in red
plot(LPa(:),'ro')
plot(LMa(:),'bo')
%plot(asymP,'go') % Asymptotic guess in green

hold off

axis([xp(1), xp(end), yp(end), yp(1)])
axis([-10,10,-10,10])
caxis([-100,-90])
axis on
%axis equal
shading interp
colorbar
drawnow
end

KMGM0=Kminus(GM0,Kargs);
KPGM0=Kplus(GM0,Kargs);
KMPM0=Kminus(PM0,Kargs);
KPPM0=Kplus(PM0,Kargs);

%% Define K derivatives
[output.KpprTM,output.KmprTP]=computeKDerivatives(ADData,AAData,TMd,KMTM,TPd,KPTP);

%% Put data into new structure
output.LMa =    LMa;
output.LPa =    LPa;
output.ZMa =    ZMa;
output.ZPa =    ZPa;
output.Kargs=   Kargs;
output.TMd =    TMd;
output.TPd =    TPd;
output.zGM0 = mysqrt(omega*w,GM0);
output.zPM0 = mysqrt(omega*w,PM0);
output.KMTM=    KMTM;
output.KPTP=    KPTP;
output.KMLM=    KMLM;
output.KPLP=    KPLP;
output.KMPM0=    KMPM0;
output.KPPM0=    KPPM0;
output.KMGM0=    KMGM0;
output.KPGM0=    KPGM0;
output.spac=    ADData.spac;
output.Beta=    ADData.Beta;
output.Sigma=   AAData.Sigma;
output.w =      AAData.w;
output.SQRTa=   SQRTa;
output.comb = Modes.comb;
end

function [KpprTM,KmprTP]=computeKDerivatives(ADData,AAData,TMd,KMTM,TPd,KPTP)
% Compute the derivatives of K.  
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
KprTM=  (TMnum./TMden.^2 + 1i*mu(2) + 2*mu(3)*TMd)/4/pi;

KpprTM=KprTM./KMTM;

TPden = cos(s*zTP)-cos(d*TPd+Sigma);
TPnum = (-TPd./zTP.*sin(s*zTP)-s*TPd.*cos(s*zTP)).*TPden...
        -zTP.*sin(s.*zTP).*(TPd.*s./zTP.*sin(s*zTP)+d*sin(d*TPd+Sigma));
KprTP=  (TPnum./TPden.^2 + 1i*mu(2) + 2*mu(3)*TPd )/4/pi;

KmprTP=KprTP./KPTP;

end