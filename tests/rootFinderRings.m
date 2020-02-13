%% This finds the roots through the ring algorithm.
addpath('cfe')

omega5 = 30; he = 1; w = 1; sigma = 3*pi/4; d =1;
mu = [20 0 0]; s = sqrt(d^2 + he^2);

AAData.omega5 = omega5;
AAData.w = w;
ADData.spac(1) = he;
ADData.spac(2) = d;
AAData.sigma5 = sigma;
ADData.mu = mu;

zeta = @(xVar) mysqrt(omega5*w,xVar);
den = @(xVar)  cos(he*zeta(xVar))-cos(d*xVar+sigma);
num = @(xVar) zeta(xVar).*sin(he*zeta(xVar)) + (mu(1) + xVar*mu(2) + xVar.^2*mu(3)).*den(xVar);

tol = 1e-10;

Modes.trunc = 1e2;
profile on

f = @(xVar) 1./KNumlogD(xVar,ADData,AAData);
logD = @(xVar) KNumlogD(xVar,ADData,AAData);

asympGuess = computeAsympGuess(ADData,AAData,Modes);
knownRootsInit = newtonKRoots(asympGuess,f,logD,tol,[]);
knownRootsInit = myUnique(knownRootsInit,tol);

nTol = 5; % The maximum number of zeros permitted in each annulus

R0 = max(2*w*omega5,max(pi/s*abs(mu(1))));
%knownRootsInit = knownRootsInit(abs(knownRootsInit)<1.1*R0);

chi = 1.3;
if cos(chi)<0; error('The ellipse does not have the right parameters'); end

tic
finRoots = rootFinder(f,logD,knownRootsInit,R0,chi,tol);
toc

%% Plots

