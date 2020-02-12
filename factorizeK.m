function factorizeK(ADData,AAData,Modes)

s = ADData.spac(1);
d = ADData.spac(2);
k = AAData.omega;
w = AAData.w;
mu = ADData.mu;
sigma = AAData.sigma;

zeta = @(xVar) mysqrt(k*w,xVar);

fn  = @(x) (zeta(x).*sin(s*zeta(x))+ mu*(cos(s*zeta(x)) - cos(d*x + sigma)));
fnp = @(x) -x./zeta(x).*sin(s*zeta(x))-s*x.*cos(s*zeta(x)) + mu*(s*x./zeta(x).*sin(s*zeta(x)) + d*sin(d*x+sigma));
argFun = @(xVar) fnp(xVar)./fn(xVar);

boxLength=30+1;
boxHeight=60;

ways = [(boxLength-1i*boxHeight)/2,(boxLength+1i*boxHeight)/2,(-boxLength+1i*boxHeight)/2];
Npr = (real(integral(argFun,-(boxLength+1i*boxHeight)/2,-(boxLength+1i*boxHeight)/2,'Waypoints',ways,'ArrayValued',true)/(2i*pi)));
N = round(Npr);
N3 = permute(1:N,[1,3,2]);

polyFun = @(xVar) bsxfun(@power,xVar,N3);
newFun = @(zVar) polyFun(zVar).*argFun(zVar);
tic
p = integral(newFun,-(boxLength+1i*boxHeight)/2,-(boxLength+1i*boxHeight)/2,'Waypoints',ways,'ArrayValued',true)/(2i*pi);

e = ones(1,N+1);

for l = 1 :N
    j = (1:l).';
    e(l+1) = 1./l.*sum(permute((-1).^(j-1).'.*e(l-j+1),[1,3,2]).*p(:,:,(1:l)),3);
end

myCloseRoots = roots((-1).^(1:N+1).*e);

closeRootsPlus = myCloseRoots(imag(myCloseRoots)>0).';
closeRootsMinus = myCloseRoots(imag(myCloseRoots)<0).';

a  = (d+1i*s)*2*pi/(d^2+s^2);
b  = 1i*(d+1i*s)/(d^2+s^2);
c  = (pi-sigma+1i*log(2*pi/mu/sqrt(d^2+s^2))-atan(s/d))*(d+1i*s)./(d^2+s^2);

a1 = (-d+1i*s)*2*pi/(d^2+s^2);
b1 =  -1i*(-d+1i*s)/(d^2+s^2);
c1 = (-pi+sigma-1i*log(2*pi/mu/sqrt(d^2+s^2))-atan(s/d))*(-d+1i*s)./(d^2+s^2);

trunc = Modes.trunc;
n = 1:trunc;

Rsols1b = a*n+b*log(n)+c;
Rsols2b = a1*n+b1*log(n)+c1;
Rsols3b = conj(a1)*n+conj(b1)*log(n)+conj(c1);
Rsols4b = conj(a)*n+conj(b)*log(n)+conj(c);

%  sols1ab = findRoots2(s,d,sigma,mu,k,w,Rsols1b);
%  sols2ab = findRoots2(s,d,sigma,mu,k,w,Rsols2b);
%  sols3ab = conj(sols2ab);
%  sols4ab = conj(sols1ab);
sols1ab = [];
sols2ab = [];
sols3ab = [];
sols4ab = [];

for nL = 1:numel(Rsols1b)
%  sols1ab = [sols1ab,findRoots(s,d,sigma,mu,k,w,[real(Rsols1b(nL));imag(Rsols1b(nL))])];
%  sols2ab = [sols2ab,findRoots(s,d,sigma,mu,k,w,[real(Rsols2b(nL));imag(Rsols2b(nL))])];
%  sols3ab = [sols3ab,findRoots(s,d,sigma,mu,k,w,[real(Rsols3b(nL));imag(Rsols3b(nL))])];
%  sols4ab = [sols4ab,findRoots(s,d,sigma,mu,k,w,[real(Rsols4b(nL));imag(Rsols4b(nL))])];
 sols1ab = [sols1ab,findRoots2(s,d,sigma,mu,k,w,Rsols1b(nL))];
 sols2ab = [sols2ab,findRoots2(s,d,sigma,mu,k,w,Rsols2b(nL))];
 sols3ab = [sols3ab,findRoots2(s,d,sigma,mu,k,w,Rsols3b(nL))];
 sols4ab = [sols4ab,findRoots2(s,d,sigma,mu,k,w,Rsols4b(nL))];
end

% sols1ab = complex(sols1ab(1,:),sols1ab(2,:));
% sols2ab = complex(sols2ab(1,:),sols2ab(2,:));
% sols3ab = complex(sols3ab(1,:),sols3ab(2,:));
% sols4ab = complex(sols4ab(1,:),sols4ab(2,:));
sols1ab(imag(sols1ab)<boxHeight/2)=[];
sols2ab(imag(sols2ab)<boxHeight/2)=[];
sols3ab(imag(sols3ab)>-boxHeight/2)=[];
sols4ab(imag(sols4ab)>-boxHeight/2)=[];

TP3 = sort(permute([closeRootsPlus,sols1ab,sols2ab],[1,3,2]));
TM3 = sort(permute([closeRootsMinus,sols3ab,sols4ab],[1,3,2]));

% Make entries unique
[~,indPr] = unique(real(TP3));
[~,indPi] = unique(imag(TP3));
indCompP = intersect(indPr,indPi);
TP3 = TP3(indCompP);

[~,indMr] = unique(real(TM3));
[~,indMi] = unique(imag(TM3));
indCompM = intersect(indMr,indMi);
TM3 = TM3(indCompM);


% plot(permute(TP3,[3,2,1]),'.')
% hold on
% plot(Rsols3b,'.')
%  hold off
%return
% plot(abs(permute(K(TM3,ADData,AAData),[3,2,1])),'.')
%plot(abs(K(closeRootsPlus,ADData,AAData)),'.')

%plot(abs(sols2ab),abs(K(sols2ab,ADData,AAData)),'.')
%return

%% Extract data from structs
%h=ADData.spac(1); d=ADData.spac(2); chie=ADData.chie;

chie=atan(d/s);
se = sqrt(d^2+s^2);
nx = 200; ny = 200;
xp = linspace(-15,15,nx);
yp = linspace(-2,40,ny);
[X,Y]=meshgrid(xp,yp);
Z = X + 1i*Y;
%gammav = linspace(2,10)*exp(1i*pi/4);
gammav = Z;
%% Define acoustic modes
aTrunc = permute(-trunc:trunc,[1,3,2]);
f=(bsxfun(@minus,sigma,2*pi*aTrunc))/se;
SQRT=mysqrt(w*k,f); %output.SQRT=SQRT;
LM3=-f*sin(chie)-cos(chie)*SQRT;
LP3=-f*sin(chie)+cos(chie)*SQRT;

%% Define numerator
numM=1-bsxfun(@rdivide,gammav,TP3);
numprodM=prod(numM,3);
%% Define denominator
denM=1-bsxfun(@rdivide,gammav,LP3);
denprodM=prod(denM,3);
%% Define constant and exponential term
Kminus=bsxfun(@times,exp((1i*gammav/pi)*(s*log(2*cos(chie))+chie*d)),numprodM./denprodM);

%% Define numerator
numP=1-bsxfun(@rdivide,gammav,TM3);
numprodP=prod(numP,3);

%% Define denominator
denP=1-bsxfun(@rdivide,gammav,LM3);
denprodP=prod(denP,3);

%% Define constant and exponential term
const = k*w.*sin(k*w*s)./(4*pi*(cos(k*w*s)-cos(sigma)))+mu/(4*pi);
%Jplus=bsxfun(@times,const,exp((-1i*gammav/pi)*(h*log(2*cos(chie))+chie*d))).*numprod./denprod;
Kplus=const.*bsxfun(@times,exp((-1i*gammav/pi)*(s*log(2*cos(chie))+chie*d)),numprodP./denprodP);

% plot(permute(cos(h*mysqrt(k*w,LP3(1:10))),[3,2,1]))
% hold on
% plot(permute(cos(d*LP3(1:10) + sigma),[3,2,1]))
% hold off
%(cos(s*mysqrt(k*w,LP3(1:2))) - cos(d*LP3(1:2) + sigma))
%K(LP3(1:2),ADData,AAData)

%plot(permute(K(TP3,ADData,AAData),[3,1,2]))

%plot(permute(K(TM3,ADData,AAData),[3,1,2]))
%plot(abs(permute(K(LM3(1:10),ADData,AAData),[3,2,1])))



figure(1)
clf
%PhasePlot(Z,K(Z,ADData,AAData)./(Kminus.*Kplus),'d');
%PhasePlot(Z,K(Z,ADData,AAData)./(numprodM./denprodM.*numprodP./denprodP),'d');
pcolor(real(Z),imag(Z),abs(K(Z,ADData,AAData)./(numprodM./denprodM.*numprodP./denprodP.*const)));
 colorbar
 caxis([.9,1.1])
shading interp;

hold on
allWays = [-(boxLength+1i*boxHeight)/2,ways,-(boxLength+1i*boxHeight)/2];
plot(permute(TP3,[3,2,1]),'wx')
plot(permute(TM3,[3,2,1]),'wx')
plot(permute(LP3,[3,2,1]),'bo')
%plot(permute(TM3,[3,2,1]),'wx')
plot(Rsols1b,'r.')
plot(allWays,'k','LineWidth',3)
axis([xp(1),xp(end),yp(1),yp(end)])
%zlim([-10,10])
%caxis([-10,10])
axis on
hold off


% figure(2)
% %PhasePlot(Z,K(Z,ADData,AAData)./(Kminus.*Kplus),'d');
% PhasePlot(Z,K(Z,ADData,AAData)./(1./denprodM),'d');
% %surf(real(Z),imag(Z),real(K(Z,ADData,AAData).*(denprodP.*denprodM)));
% shading interp;
% 
% hold on
% allWays = [-(boxLength+1i*boxHeight)/2,ways,-(boxLength+1i*boxHeight)/2];
% plot(permute(TP3,[3,2,1]),'wx')
% plot(permute(TM3,[3,2,1]),'wx')
% plot(permute(LP3,[3,2,1]),'bo')
% %plot(permute(TM3,[3,2,1]),'wx')
% plot(Rsols1b,'r.')
% plot(allWays,'k','LineWidth',3)
% axis([xp(1),xp(end),yp(1),yp(end)])
% zlim([-10,10])
% caxis([-10,10])
% axis on
% hold off
%return
%Kplus.*Kminus - K(0,ADData,AAData);
% %return
% plot(abs(gammav),Kplus.*Kminus)
% hold on
% plot(abs(gammav),K(gammav,ADData,AAData))
% hold off
% %return
% Kplus.*Kminus./K(gammav,ADData,AAData);
%return
figure(2)
plot(abs(gammav),abs(Kplus.*Kminus./K(gammav,ADData,AAData)))
hold on
%plot(abs(gammav),K(gammav,ADData,AAData))
ylim([0,2])
hold off

figure(3)

plot(imag(permute(TP3,[3,2,1])),'o')
return
figure(2)
K(sols1ab,ADData,AAData);
plot(abs(K(closeRootsPlus,ADData,AAData)))

end