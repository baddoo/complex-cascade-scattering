
function phi = computeField(Z,data,type)
%COMPUTEFIELD Computes that field at the points Z

spac = data.spac;
s = spac(1); d = spac(2); del = sqrt(s^2 + d^2);
LPa = data.LPa(:).'; % Row vector of Lambda^+
ZPa = data.ZPa(:).'; % Row vector of Zeta^+
LMa = data.LMa(:).'; % Row vector of Lambda^+
ZMa = data.ZMa(:).'; % Row vector of Zeta^+
SQRT = data.SQRTa(:).'; % Row vector of SQRT    
ADData = data.ADData;
AAData = data.AAData;
TMd = data.TMd(:).';
TPd = data.TPd(:).';
GM = data.GM0;
PM = data.PM0;
omega = data.AAData.omega;
w = data.w;
zeta = @(zVar) mysqrt(omega*w,zVar);
A = data.A(:);
B = data.B(:);
C = data.C(:);
KPGM = data.KPGM0;
KMPM = data.KMPM0;

zetaGM0 = zeta(GM);
zetaPM0 = zeta(PM);
w0 = data.AAData.Amp(2);

M = data.ADData.M;
Sigma = data.Sigma;

% Partition Z into different sections 
Zsz = size(Z); % Original size of Z
Z = Z(:); % Turn Z into a column vector
X = real(Z); Y = imag(Z);

% Calculate upstream field
data.comb=[1,0,1,0];
upLoc = (X <= Y*d/s);
Xup = X(upLoc); Yup = Y(upLoc);  % Upsteam region
D13LP = D(LPa,data); %  D^{(1,3)} evaluated at Lambda^+
ARLP = AdR(Xup,Yup,LPa,ZPa,SQRT,spac,PM,M,type) + AuR(Xup,Yup,LPa,ZPa,SQRT,spac,PM,M,type); % 
phiUp = 2i*pi*ARLP*D13LP.';
% Alternative expression of potential solution
%phiUp = pi/del*(ZPa.*exp(1i*ZPa.*Yup).*exp(-1i*LPa.*Xup))*(D13LP./SQRT).';

% Calculate upper triangular field
utLoc = (X > Y*d/s & X<=d);
Xut = X(utLoc); Yut = Y(utLoc); % Upstream triangle

phiUt =  2i*pi*Ad(Xut,Yut,TPd,zeta(TPd),spac,Sigma,PM,M,type)*B ...
       - 2*pi*Ad(Xut,Yut,GM,zetaGM0,spac,Sigma,PM,M,type)*w0/(2*pi)^2/regK(GM,ADData,AAData) ...
       - 2*pi*Ad(Xut,Yut,TMd,zeta(TMd),spac,Sigma,PM,M,type)./(TMd - PM)*(A + C) ...
       + 2i*pi*AuR(Xut,Yut,LPa,ZPa,SQRT,spac,PM,M,type)*D13LP.';

% Calculate rectangular field
reLoc = (X>d & X<=2);
Xre = X(reLoc); Yre = Y(reLoc);  % Middle rectangle
ATMd = Ad(Xre,Yre,TMd,zeta(TMd),spac,Sigma,PM,M,type) + Au(Xre,Yre,TMd,zeta(TMd),spac,Sigma,PM,M,type);
ATPd = Ad(Xre,Yre,TPd,zeta(TPd),spac,Sigma,PM,M,type) + Au(Xre,Yre,TPd,zeta(TPd),spac,Sigma,PM,M,type);
AGM0 = Ad(Xre,Yre,GM,zeta(GM),spac,Sigma,PM,M,type) + Au(Xre,Yre,GM,zeta(GM),spac,Sigma,PM,M,type);

phiRe =  2i*pi*ATPd*B ...
       - AGM0*w0/(2*pi)/regK(GM,ADData,AAData) ...
       - 2*pi*ATMd./(TMd - PM)*(A + C);

% Calculate downstream triangular field
% Calculate wake coefficient
AsumG = -sum((A+C)./(1i*(TMd.'-PM)).*data.KMTM(:)./data.KMPM0.*exp(-2i*(TMd.'-PM)));
TtermsG = -w0/(4i*pi^2)./KPGM./KMPM.*exp(-2i*(GM-PM));
P =  AsumG + TtermsG;
dtLoc = ( X <2+Y*d/s & X > 2);
Xdt = X(dtLoc); Ydt = Y(dtLoc); % Downstream triangle
data.comb=[0,1,0,1];
D24LM = D(LMa,data); %  D^{(1,3)} evaluated at Lambda^+
phiDt =  2i*pi*Au(Xdt,Ydt,TPd,zeta(TPd),spac,Sigma,PM,M,type)*B ...
       - 2*pi*Au(Xdt,Ydt,GM,zetaGM0,spac,Sigma,PM,M,type)*w0/(2*pi)^2/regK(GM,ADData,AAData) ...
       - 2*pi*Au(Xdt,Ydt,TMd,zeta(TMd),spac,Sigma,PM,M,type)./(TMd - PM)*(A + C) ...
       - 2i*pi*AdR(Xdt,Ydt,LMa,ZMa,SQRT,spac,PM,M,type)*D24LM.' ...
       + 2i*pi*P*Ad(Xdt,Ydt,PM,zetaPM0,spac,Sigma,PM,M,type);

% Calculate downstream field
dsLoc = ( X > 2 + Y*d/s );
Xdn = X(dsLoc); Ydn = Y(dsLoc);  % Downstream region
ARLM = AdR(Xdn,Ydn,LMa,ZMa,SQRT,spac,PM,M,type) + AuR(Xdn,Ydn,LMa,ZMa,SQRT,spac,PM,M,type); 
APM = Ad(Xdn,Ydn,PM,zeta(PM),spac,Sigma,PM,M,type) + Au(Xdn,Ydn,PM,zeta(PM),spac,Sigma,PM,M,type);
phiDs = -2i*pi*ARLM*D24LM.' ...
        + 2i*pi*P*APM;
% Alternative expression of potential solution
% phiDs = - pi/del*(ZMa.*exp(-1i*ZMa.*Ydn).*exp(-1i*LMa.*Xdn))*(D24LM./SQRT).' ...
%         + 2i*pi*P*APM;

% Now put components back into the correct entries in the matrix
phi = nan(Zsz);
phi(upLoc) = phiUp;
phi(utLoc) = phiUt;
phi(reLoc) = phiRe;
phi(dtLoc) = phiDt;
phi(dsLoc) = phiDs;

end