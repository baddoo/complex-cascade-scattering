mu = [20,2,10i--1i];
ADData.mu = mu;

he = 1;
d = .25;

ADData.spac = [he,d];
s = sqrt(he^2+d^2);
sigma5 = 1;
AAData.sigma5 = sigma5;
w = 1;
k5 = 10+0.0i;
AAData.w = w;
AAData.k5 = k5;
trunc = 10;

Modes.trunc = trunc;
asympGuess = computeAsympGuess(ADData,AAData,Modes);

asympRoots = [];

for nL = 1:numel(asympGuess)
    asympRoots = [asympRoots,findRoots3(he,d,sigma5,mu,k5,w,asympGuess(nL))];
    %sols2ab = [sols2ab,findRoots3(he,d,sigma5,mu,k5,w,Rsols2b(nL))];
    %sols3ab = [sols3ab,findRoots3(he,d,sigma5,mu,k5,w,Rsols3b(nL))];
    %sols4ab = [sols4ab,findRoots3(he,d,sigma5,mu,k5,w,Rsols4b(nL))];
end

figure(1)
loglog(abs(asympGuess-asympRoots).')
%plot(abs(Rsols1b))
hold on
%plot(abs(sols1ab))
%loglog(abs(Rsols2b-sols2ab))
% loglog(abs(Rsols3b-sols3ab))
% loglog(abs(Rsols4b-sols4ab))


hold off

nx = 300; ny = 300;
xp = linspace(-20,20,nx);
yp = -linspace(-70,70,ny);
[X,Y]=meshgrid(xp,yp);
Z = X + 1i*Y;
zeta = @(xVar) mysqrt(k5*w,xVar);

fn  = @(x) zeta(x).*sin(he*zeta(x)) + ...
           (mu(1)+mu(2).*x + mu(3).*x.^2).*(cos(he*zeta(x)) - cos(d*x + sigma5));
fd = @(x) cos(he*zeta(x)) - cos(d*x + sigma5);
max(abs(fn(asympRoots)./fd(asympRoots)))
pl=0;
figure(2)
clf
PhasePlot(Z,fn(Z),'d');
% 
hold on
plot(asympGuess,'gx')
plot(asympRoots,'kx')
plot([asympGuess;asympRoots])
%plot(asymP,'go')

hold off

axis([xp(1), xp(end), yp(end), yp(1)])
axis on
shading interp
colorbar
