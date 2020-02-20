% Finds the roots of the meromorphic function.


function finRoots = rootFinder(f,logD,knownRootsInit,R0,chi,tol)

nTol = 5;

polyLogD = @(xVar,krV) sum(1./(xVar - permute(krV(:),[4,2,3,1])),4);
argFun = @(xVar,Np,krV) (logD(xVar) - polyLogD(xVar,krV)).*xVar.^Np;

% Define ellipse
z = @(th,rV) rV.*(cos(th) + 1i*cos(chi)*sin(th));
dzdth = @(th,rV) rV.*(-sin(th) + 1i*cos(chi)*cos(th));

% Define integrand
intFun = @(th,rV,Np,krV) dzdth(th,rV).*argFun(z(th,rV),Np,krV);

% Define roots that are already known
newKR = knownRootsInit;

% Count zeros in R0-circle
nInt = 1e2*round(2*pi*R0);
tInt = linspace(0,2*pi,nInt+1); tInt(end) = [];
tInt0 = tInt; nInt0 = nInt;
integrandError = norm(intFun(tInt0,R0,0,newKR),'inf');
counter = 1;
while integrandError*abs(tInt0(1)-tInt0(2))>.1 || isnan(integrandError)
       R0 = R0+.05*counter; % INCREASE R
       nInt0 = round((1.05.^counter)*nInt0); % Increase number of quadrature nodes exponentially 
       tInt0 = linspace(0,2*pi,nInt0+1); tInt0(end) = []; % Calculate new quadrature points
       integrandError = norm(intFun(tInt0,R0,0,newKR),'inf'); 
       counter = counter + 1; % Increase counter by 1
       if counter>15; warning('The first root-counting integral is taking a long time to compute. Consider using different parameters.'); end
end

err = inf;
m = 0;
while err~=0
Rvec = [R0,0];    
NR = round(-1i/nInt0*sum(intFun(tInt0,R0,0,newKR),2));  
aNR = NR;
tInt = linspace(0,2*pi,nInt+1); tInt(end) = [];

while any(aNR > nTol)

% Find annulus with too many roots
loc = find(aNR>nTol,1);

% Find average radius
newR = (Rvec(loc) + Rvec(loc+1))/2;

% Compute number of roots in new circle

% If the function is too big on the integration contour, move the radius
% inwards whilst simultaneously increasing the number of quadrature points.
% Also, if the number of zeros in the new disc is lower than the number of zeros 
% in the next biggest disc, we know that there has been an error. In that
% case, increase the number of quadrature points until the number of zeros
% identified is smaller than that in a larger disc.
counter = 1;
intVal = intFun(tInt,newR,0,newKR);
newN = round(-1i/nInt*sum(intVal,2));
integrandError = norm(intVal,'inf');

while integrandError*abs(tInt(1)-tInt(2))>1 || isnan(integrandError) || newN>NR(loc)
   newR = newR-.1*(newR-Rvec(loc+1)); % Move the new radius closer to the interior radius 
   nInt = round(1.05.^counter*nInt); % Increase number of quadrature points by 1 percent 
   tInt = linspace(0,2*pi,nInt+1); tInt(end) = [];
   intVal = intFun(tInt,newR,0,newKR);
   integrandError = norm(intVal,'inf');
   newN = round(-1i/nInt*sum(intVal,2));
   counter = counter + 1;
end

% Place new # of roots and radius in vectors
NR = [NR(1:loc), newN, NR(loc+1:end)];
Rvec = [Rvec(1:loc), newR, Rvec(loc+1:end)];

% Calculate number of zeros in each annulus
aNR = [-diff((NR)),NR(end)];
end

% Find annuli where the there are no roots and remove outer radius from
% vectors
zerLoc = find(aNR == 0);
aNR(zerLoc) = [];
NR(zerLoc) = [];
Rvec(zerLoc) = [];

% Calculate Newton sums
mat = repmat(1:max(aNR),[numel(Rvec)-1,1]);  % remove one from R because of zero radiues
mat34 = permute(mat,[3,4,1,2]);
R3 = permute(Rvec,[1,3,2]);

p = -1i/nInt*sum(intFun(tInt,R3(1:end-1),mat34,newKR),2);
pr = permute(p,[3,4,1,2]);
pd = [-diff(pr,1,1);pr(end,:)].';

% Find the number of rings in each annulus
initRoots = newtonIdentitySolver(pd,aNR);
initRoots(initRoots==0)=[];

% Hone in on roots
refinedRoots = newtonKRoots(initRoots,f,logD,tol,newKR);

% The total number of roots in the region is
err = NR(1)-numel(myUnique(refinedRoots,tol));
NR = err;
newKR = [newKR(:);myUnique(refinedRoots(:),tol)];

% Double the number of quadrature points in the next integration.
nInt = 2*nInt;
m = m + 1;

end

finRoots = newKR(:);

% % Uncomment to view plots
% 
% nx = 1e3; ny = 1e3;
% zp = 1.1*(linspace(-max(Rvec),max(Rvec),nx) + 1i*linspace(-max(Rvec),max(Rvec),ny).');
% tic
% PhasePlot(zp,f(zp));
% toc
% hold on
% plot(initRoots,'ko','MarkerFaceColor','w','MarkerSize',10)
% plot(finRoots,'kp','MarkerFaceColor','w','MarkerSize',10)
% plot(knownRootsInit,'ks','MarkerFaceColor','w','MarkerSize',10)
% plot(z(linspace(0,2*pi)',Rvec),'LineWidth',3,'Color',[0 0 0 .25]);
% hold off
% 
% axis equal
% axis([min(real(zp(:))),max(real(zp(:))),min(imag(zp(:))),max(imag(zp(:)))])

end