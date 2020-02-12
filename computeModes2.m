function output = computeModes2(ADData,AAData,Modes)

%% Get data from structures
he=ADData.spac(1); d=ADData.spac(2); se = ADData.spac(3);
trunc=Modes.trunc; w=AAData.w; sigma5=AAData.sigma5;
k5=AAData.k5; omega5 = AAData.omega5; chie = ADData.chie; delta = ADData.delta;
amodes=Modes.amodes;
%% Define duct mods
R0 = 2*max([10,w*omega5,pi/se*abs(ADData.mu(1))]);
chi = 0;
if cos(chi)<0; error('The ellipse does not have the right parameters'); end
[TP,TM,asymP] = findDuctModes2(R0,chi,ADData,AAData,Modes);

%% Define acoustic modes
aTrunc = permute(-trunc:trunc,[1,3,2]);
f=(bsxfun(@minus,sigma5,2*pi*aTrunc))/se;
SQRT=mysqrt(w*omega5,f); output.SQRT=SQRT;
LM=-f*sin(chie)-cos(chie)*SQRT;
LP=-f*sin(chie)+cos(chie)*SQRT;

%% Define
dmodes = Modes.dmodes;

TMd = TM(:,:,(1:dmodes),:,:);
TPd = TP(:,:,(1:dmodes),:,:);

ZM= -(sigma5+d*LM-2*pi*permute(-trunc:trunc,[1,3,2,4,5]))/he;
ZP=  (sigma5+d*LP-2*pi*permute(-trunc:trunc,[1,3,2,4,5]))/he;

LMa= LM(:,:,(ceil(end/2)-amodes):(ceil(end/2)+amodes),:,:); 
LPa= LP(:,:,(ceil(end/2)-amodes):(ceil(end/2)+amodes),:,:); 
ZMa= ZM(:,:,(ceil(end/2)-amodes):(ceil(end/2)+amodes),:,:); 
ZPa= ZP(:,:,(ceil(end/2)-amodes):(ceil(end/2)+amodes),:,:); 

SQRTa = SQRT(:,:,(ceil(end/2)-amodes):(ceil(end/2)+amodes),:,:);

Kargs.ADData=ADData; Kargs.AAData=AAData;
Kargs.TM=TM; Kargs.LM=LM; Kargs.TP=TP; Kargs.LP=LP;

KMTM=permute(Kminus(permute(TMd,[3,2,1,4,5]),Kargs),[3,2,1,4,5]);
KPTP=permute(Kplus(permute(TPd,[3,2,1,4,5]),Kargs),[3,2,1,4,5]); 

KMLM=permute(Kminus(permute(LMa,[3,2,1,4,5]),Kargs),[3,2,1,4,5]);
KPLP=permute(Kplus(permute(LPa,[3,2,1,4,5]),Kargs),[3,2,1,4,5]); 

GM0 = -real(k5)*delta;
PM0 = -real(omega5)*delta;

output.GM0 = GM0;
output.PM0 = PM0;

largeVal = 20*exp(1i*linspace(-pi,pi,20));
error = max(abs(Kminus(largeVal,Kargs).*Kplus(largeVal,Kargs)./regK(largeVal,ADData,AAData)-1));
disp(['The approximate error is ', num2str(error)]);

pl=0;
if pl == 1
    figure(1)

nx = 200; ny = 200;
xp = linspace(-50,50,nx);
yp = -linspace(-50,50,ny);
[X,Y]=meshgrid(xp,yp);
Z = X + 1i*Y;
PhasePlot(Z,regK(Z,ADData,AAData)./(Kplus(Z,Kargs).*Kminus(Z,Kargs)),'d');

hold on
plot(permute(TM,[3,2,1]),'gx')
plot(permute(TP,[3,2,1]),'kx')
plot(asymP,'go')

hold off

axis([xp(1), xp(end), yp(end), yp(1)])
axis on
shading interp
colorbar

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
output.zGM0 = mysqrt(omega5*w,GM0);
output.zPM0 = mysqrt(omega5*w,PM0);

output.KMTM=    KMTM;
output.KPTP=    KPTP;
output.KMLM=    KMLM;
output.KPLP=    KPLP;
output.KMPM0=    KMPM0;
output.KPPM0=    KPPM0;
output.KMGM0=    KMGM0;
output.KPGM0=    KPGM0;

output.comb=    Modes.comb;
output.dmodes=  dmodes;
output.amodes=  Modes.amodes;
output.TP =     TP;
output.TM=      TM;
output.spac=    ADData.spac;
output.Beta=    ADData.Beta;
output.sigma5=   AAData.sigma5;
output.w =      AAData.w;
output.SQRTa=   SQRTa;

end