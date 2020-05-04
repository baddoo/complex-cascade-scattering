%% Computes that field at the relevant points

function phi = computeField(data,type)

phi = @(Z) evaluateField(Z,data,type);

function phi = evaluateField(Z,data,type)

spac = data.spac;
s = spac(1); d = spac(2);
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
omega = data.omega;
w = data.w;
zeta = @(zVar) mysqrt(omega*w,zVar);
A = data.D1.A(:);
B = data.D1.B(:);
C = data.D1.C(:);

zetaGM0 = zeta(GM);
zetaPM0 = zeta(PM);
w0 = data.AAData.Amp(2);

M = data.M;
Sigma = data.Sigma;

% Partition Z into different sections 
Zsz = size(Z); % Original size of Z
Z = Z(:); % Turn Z into a column vector
X = real(Z); Y = imag(Z);

% % Define multiplication factors to allow for different fields
% if strcmp(type,'potential')
% tFacPM = 1;
% tFacGM = 1;
% tFacLM = 1;
% tFacLPu = 1;
% tFacLPt = 1;
% tFacTM = 1;
% tFacTP = 1;
% finFac = 1;
% elseif strcmp(type,'vvelocity')
% tFac = exp(1i*X*PM*M^2);
% tFacPM = 1;
% tFacGM = 1;
% tFacLM = 1;
% tFacLPu = 1;
% tFacLPu = 1i*ZPa;
% tFacLPt = 1;
% tFacTM = 1;
% tFacTP = 1;
% finFac = 1;
% elseif strcmp(type,'hvelocity')
% tFac = 1;
% elseif strcmp(type,'pressure')
% tFacGM = GM - PM;
% tFacPM = 0;
% tFacLM = LMa - PM;
% tFacLPu = LPa - PM;
% tFacLPt = LPa - PM;
% tFacTM = TMd - PM;
% tFacTP = TPd - PM;
% finFac = -1i*exp(1i*X*PM.*M^2);
% end

% Calculate upstream field
data.comb=[1,0,1,0];
upLoc = (X <= Y*d/s);
Xup = X(upLoc); Yup = Y(upLoc);  % Upsteam region
D13LP = D(LPa,data); %  D^{(1,3)} evaluated at Lambda^+
ARLP = AdR(Xup,Yup,LPa,ZPa,SQRT,spac,PM,M,type) + AuR(Xup,Yup,LPa,ZPa,SQRT,spac,PM,M,type); % 
phiUp = pi*1i*ARLP*D13LP.';

% Calculate upper triangular field
utLoc = (X > Y*d/s & X<=d);
Xut = X(utLoc); Yut = Y(utLoc); % Upstream triangle

phiUt =  pi*1i*Ad(Xut,Yut,TPd,zeta(TPd),spac,Sigma,PM,M,type)*B ...
       - pi*Ad(Xut,Yut,GM,zetaGM0,spac,Sigma,PM,M,type)*w0/(2*pi)^2/regK(GM,ADData,AAData) ...
       - pi*Ad(Xut,Yut,TMd,zeta(TMd),spac,Sigma,PM,M,type)./(TMd - PM)*(A + C) ...
       + pi*1i*AuR(Xut,Yut,LPa,ZPa,SQRT,spac,PM,M,type)*D13LP.';

% Calculate rectangular field
reLoc = (X>d & X<=2);
Xre = X(reLoc); Yre = Y(reLoc);  % Middle rectangle
ATMd = Ad(Xre,Yre,TMd,zeta(TMd),spac,Sigma,PM,M,type) + Au(Xre,Yre,TMd,zeta(TMd),spac,Sigma,PM,M,type);
ATPd = Ad(Xre,Yre,TPd,zeta(TPd),spac,Sigma,PM,M,type) + Au(Xre,Yre,TPd,zeta(TPd),spac,Sigma,PM,M,type);
AGM0 = Ad(Xre,Yre,GM,zeta(GM),spac,Sigma,PM,M,type) + Au(Xre,Yre,GM,zeta(GM),spac,Sigma,PM,M,type);

phiRe =  pi*1i*ATPd*B ...
       - pi*AGM0*w0/(2*pi)^2/regK(GM,ADData,AAData) ...
       - pi*ATMd./(TMd - PM)*(A + C);

% Calculate downstream triangular field
% Calculate wake coefficient
AsumG = -sum((A+C)./(1i*(TMd.'-PM)).*data.KMTM(:)./data.KMPM0.*exp(-2i*(TMd.'-PM)));
TtermsG = -data.D1.T.*data.KMGM0./data.KPGM0./data.KMPM0.*exp(-2i*(GM-PM));
P =  AsumG + TtermsG;
dtLoc = ( X <2+Y*d/s & X > 2);
Xdt = X(dtLoc); Ydt = Y(dtLoc); % Downstream triangle
data.comb=[0,1,0,1];
D24LM = D(LMa,data); %  D^{(1,3)} evaluated at Lambda^+
phiDt =  pi*1i*Au(Xdt,Ydt,TPd,zeta(TPd),spac,Sigma,PM,M,type)*B ...
       - pi*Au(Xdt,Ydt,GM,zetaGM0,spac,Sigma,PM,M,type)*w0/(2*pi)^2/regK(GM,ADData,AAData) ...
       - pi*Au(Xdt,Ydt,TMd,zeta(TMd),spac,Sigma,PM,M,type)./(TMd - PM)*(A + C) ...
       - pi*1i*AdR(Xdt,Ydt,LMa,ZMa,SQRT,spac,PM,M,type)*D24LM.' ...
       + 1i*pi*P*Ad(Xdt,Ydt,PM,zetaPM0,spac,Sigma,PM,M,type);

% Calculate downstream field
dsLoc = ( X > 2 + Y*d/s );
Xdn = X(dsLoc); Ydn = Y(dsLoc);  % Downstream region
ARLM = AdR(Xdn,Ydn,LMa,ZMa,SQRT,spac,PM,M,type) + AuR(Xdn,Ydn,LMa,ZMa,SQRT,spac,PM,M,type); 
APM = Ad(Xdn,Ydn,PM,zeta(PM),spac,Sigma,PM,M,type) + Au(Xdn,Ydn,PM,zeta(PM),spac,Sigma,PM,M,type);
phiDs = - pi*1i*ARLM*D24LM.' ...
        + 1i*pi*P*APM;

% Now put components back into the correct entries in the matrix
phi = nan(Zsz);
phi(upLoc) = phiUp;
phi(utLoc) = phiUt;
phi(reLoc) = phiRe;
phi(dtLoc) = phiDt;
phi(dsLoc) = phiDs;

end

end







