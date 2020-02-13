function [newADData,newAAData] = prepareData(ADData,AAData,Modes)

newADData = ADData;
newAAData = AAData;

M = ADData.M;
U = M*ADData.c0;
Beta = sqrt(1-M^2);
newADData.Beta = Beta;
delta = 1/Beta^2;
newADData.delta = delta;
physC = ADData.physC;


he = ADData.spac(1)*Beta; % Update space with PG transformation and non-dim
d = ADData.spac(2);

% if isempty(AAData.k)
%     omega = AAData.omega;
%     k = omega.*physC/(2*U);
% elseif isempty(AAData.omega)
k = AAData.k;
omega = AAData.omega;
%     omega=k./ADData.physC.*(2*U);    
% end

if isempty(AAData.kn)
    sigmao = AAData.sigmao;
    sigma= d*omega*delta*M^2 + sigmao;
    kn =  (sigma - k*delta*d)/(omega*he);
elseif isempty(AAData.sigmao)
    kn = AAData.kn;
    sigma = k*delta*d + omega*kn*he;
    sigmao= sigma - d*omega*delta*M^2;
end    

newAAData.kn = kn;
k = AAData.k;
omega = AAData.omega;

newADData.spac(1) = he;
newADData.spac(2) = d;

newADData.spac(3) = sqrt(d^2+he^2); % Do PG transformation
newADData.chie = atan(d/he); % Do PG transformation

newAAData.sigma5 = reshape(sigma,[1 1 1 1 numel(sigma)]);
newAAData.sigmao5 = reshape(sigmao,[1 1 1 1 numel(sigmao)]);

newAAData.k5 = reshape(k,[1 1 1 1 numel(k)]);
newAAData.omega5 = reshape(omega,[1 1 1 1 numel(omega)]);

newAAData.Amp5 = [0,1i*newAAData.omega5*kn,0];

w=mysqrt(M*delta,newAAData.kz/Beta); 
newAAData.w=w;

C = ADData.C;
if strcmp(ADData.case,'I')
    
   mu = C*[1,0,0];
   
elseif strcmp(ADData.case,'II')
    
    mu0 = 1i*omega*delta;
    mu1 = -1;
    mu = C*[mu0,mu1,0];
    
elseif strcmp(ADData.case,'III')   
    
    mu0 = -2*omega^2;
    mu1 = -2i*M^2*omega; 
    mu2 = 1;
    mu = C*[mu0,mu1,mu2];
    
end

newADData.mu = mu;

end