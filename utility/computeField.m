%% Computes that field at the relevant points

function phi = computeField(data,type)

% Extract data from structure
% 
% Define complexified coordinate
% Z = X + 1i*Y;
% Zup = Z(X>=Y*d/s); % Upsteam region
% Zut = Z(X<Y*ds/s | X>=d); % Upstream triangle
% Zre = Z(X<d | X>=2); % Middle rectangle
% Zdt = Z(X>=2 + Y*d/s | X<2); % Downstream triangle
% Zdn = Z(X>=2); % Downstream region
% 
% upstream = field(data,type);
% upperTriangle = acousticFieldUpperTriangle(data,type);
% rectangle = acousticFieldRectangle(data,type);
% lowerTriangle = acousticFieldLowerTriangle(data,type);
% downstream = acousticFieldDownstream(data,type);
% 
% m = [1,1,1,1,size(kx,5)];
% 
% 

phi = @(Z) evaluateField(Z,data);

function phi = evaluateField(Z,data)

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
GM0 = data.GM0;
PM0 = data.PM0;
omega = data.omega;
w = data.w;
zeta = @(zVar) mysqrt(omega*w,zVar);

A = data.D1.A(:);
B = data.D1.B(:);
C = data.D1.C(:);

zetaGM0 = zeta(GM0);
zetaPM0 = zeta(PM0);
w0 = data.AAData.Amp(2);

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
ARLP = AdR(Xup,Yup,LPa,ZPa,SQRT,spac) + AuR(Xup,Yup,LPa,ZPa,SQRT,spac); % 
phiUp = pi*1i*ARLP*D13LP.';

% Calculate upper triangular field
utLoc = (X > Y*d/s & X<=d);
Xut = X(utLoc); Yut = Y(utLoc); % Upstream triangle
phiUt =  pi*1i*Ad(Xut,Yut,TPd,zeta(TPd),spac,Sigma)*B ...
       - pi*Ad(Xut,Yut,GM0,zetaGM0,spac,Sigma)*w0/(2*pi)^2/K(GM0,ADData,AAData) ...
       - pi*Ad(Xut,Yut,TMd,zeta(TMd),spac,Sigma)./(TMd - PM0)*(A + C) ...
       + pi*1i*AuR(Xut,Yut,LPa,ZPa,SQRT,spac)*D13LP.';

% Calculate rectangular field
reLoc = (X>d & X<=2);
Xre = X(reLoc); Yre = Y(reLoc);  % Middle rectangle
ATMd = Ad(Xre,Yre,TMd,zeta(TMd),spac,Sigma) + Au(Xre,Yre,TMd,zeta(TMd),spac,Sigma);
ATPd = Ad(Xre,Yre,TPd,zeta(TPd),spac,Sigma) + Au(Xre,Yre,TPd,zeta(TPd),spac,Sigma);
AGM0 = Ad(Xre,Yre,GM0,zeta(GM0),spac,Sigma) + Au(Xre,Yre,GM0,zeta(GM0),spac,Sigma);

phiRe =  pi*1i*ATPd*B ...
       - pi*AGM0*w0/(2*pi)^2/K(GM0,ADData,AAData) ...
       - pi*ATMd./(TMd - PM0)*(A + C);

AsumG = -sum((A+C)./(1i*(TMd.'-PM0)).*data.KMTM(:)./data.KMPM0.*exp(-2i*(TMd.'-PM0)));
TtermsG = -data.D1.T.*data.KMGM0./data.KPGM0./data.KMPM0.*exp(-2i*(GM0-PM0));
P =  AsumG + TtermsG;

% Calculate downstream triangular field
dtLoc = ( X <2+Y*d/s & X > 2);
Xdt = X(dtLoc); Ydt = Y(dtLoc); % Downstream triangle
data.comb=[0,1,0,1];
D24LM = D(LMa,data); %  D^{(1,3)} evaluated at Lambda^+
phiDt =  pi*1i*Au(Xdt,Ydt,TPd,zeta(TPd),spac,Sigma)*B ...
       - pi*Au(Xdt,Ydt,GM0,zetaGM0,spac,Sigma)*w0/(2*pi)^2/K(GM0,ADData,AAData) ...
       - pi*Au(Xdt,Ydt,TMd,zeta(TMd),spac,Sigma)./(TMd - PM0)*(A + C) ...
       - pi*1i*AdR(Xdt,Ydt,LMa,ZMa,SQRT,spac)*D24LM.' ...
       + 1i*pi*P*Ad(Xdt,Ydt,PM0,zetaPM0,spac,Sigma);

% Calculate downstream field
dsLoc = ( X > 2 + Y*d/s );
Xdn = X(dsLoc); Ydn = Y(dsLoc);  % Downstream region
ARLM = AdR(Xdn,Ydn,LMa,ZMa,SQRT,spac) + AuR(Xdn,Ydn,LMa,ZMa,SQRT,spac); 
APM0 = Ad(Xdn,Ydn,PM0,zeta(PM0),spac,Sigma) + Au(Xdn,Ydn,PM0,zeta(PM0),spac,Sigma);
phiDs = - pi*1i*ARLM*D24LM.' ...
        + 1i*pi*P*APM0;


phi = nan(Zsz);
phi(upLoc) = phiUp;
phi(utLoc) = phiUt;
phi(reLoc) = phiRe;
phi(dtLoc) = phiDt;
phi(dsLoc) = phiDs;

end






end

% upstream(repmat(X>=Y*d/s,m))=0;
% upperTriangle(repmat(X<Y*d/s | X>=d,m))=0;
% rectangle(repmat(X<d | X>=2,m))=0;
% lowerTriangle(repmat(X>=2+Y*d/s | X<2,m))=0;
% downstream(repmat(X-2<Y*d/s,m))=0;
%









