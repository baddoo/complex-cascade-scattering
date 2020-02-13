function [newADData,newAAData] = prepareData(ADData,AAData)

newADData = ADData;
newAAData = AAData;

M = ADData.M;
Udim = M*ADData.c0;
Wdim = ADData.Wdim;
Beta = sqrt(1-M^2);
newADData.Beta = Beta;
chordDim = ADData.chordDim;
semiChordDim = chordDim/2;
omegaDim = AAData.omegaDim;
kxDim = AAData.kxDim;
kyDim = AAData.kyDim;
kzDim = AAData.kzDim;

s = ADData.spacDim(1)*Beta/semiChordDim; % Update space with PG transformation and non-dim
d = ADData.spacDim(2)/semiChordDim;

newADData.spac(1) = s;
newADData.spac(2) = d;

newADData.spac(3) = sqrt(d^2+s^2); % Do PG transformation
newADData.chie = atan(d/s); % Do PG transformation

% if isempty(AAData.k)
%     omega = AAData.omega;
%     k = omega.*physC/(2*U);
% elseif isempty(AAData.omega)
%     omega=k./ADData.physC.*(2*U);    
% end
omega = semiChordDim/Udim*(omegaDim - Wdim*kzDim);

if strcmp(kxDim,'gust')
    kx = omega;
else
    kx = kxDim*Beta^2*semiChordDim + M^2*omega;
end

ky = kyDim*semiChordDim/omega/Beta;
kz = kzDim*semiChordDim/omega;

if isempty(AAData.kyDim)
    sigmao = AAData.sigmao;
    sigma= d*omega*(M/Beta)^2 + sigmao;
    ky =  (sigma - kx/Beta^2*d)/(omega*s);
elseif isempty(AAData.sigmao)
    ky = AAData.kn;
    sigma = kx/Beta^2*d + omega*ky*s;
    sigmao= sigma - d*omega*(M/Beta)^2;
end    

newAAData.sigma = sigma;
newAAData.sigmao = sigmao;
newAAData.omega = omega;
newAAData.kx = kx;
newAAData.ky = ky;
newAAData.kz = kz;

newAAData.Amp = [0,1i*omega*ky,0];

w=mysqrt(M/Beta^2,kz/Beta); % Need to spanwise flow terms.
newAAData.w=w;

C = ADData.C;

switch ADData.case
    case 1
    
        mu = C*[1,0,0];
   
    case 2
    
        mu0 = 1i*omega/Beta^2;
        mu1 = -1;
        mu = C*[mu0,mu1,0];
    
    case 3
    
        % Need to complete    
        mu0 = -2*omega;
        mu1 = -2i*M^2*omega; 
        mu2 = 1;
        mu = C*[mu0,mu1,mu2];
    
end

newADData.mu = mu;

end