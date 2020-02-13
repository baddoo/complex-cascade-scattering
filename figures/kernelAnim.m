omega = 5*(1+1e-3i); w = 1; dPG = .5; sigma = 1;
mu = [1i,0,0];
addpath('cfe')
x = linspace(-15,15,1e3);
y = linspace(-15,15,1e3);
z = x+ 1i*y.';


figure(1)
PhasePlot(z,mysqrt(omega*w,z),'e'); % or e
h = gca;
set(h,'position',[0 0 1 1]);
print('../../../general-presentations/images/animated-images/leppington','-dpng','-r50')
%%
sVec = 2.^linspace(0,5,180);
for l = 1:numel(sVec)
sPG = sVec(l);
F= @(x) -1i*mysqrt(omega*w,x).*sin(sPG*mysqrt(omega*w,x))./(cos(sPG*mysqrt(omega*w,x)) - cos(dPG*x + sPG*sigma)) ...
        + (-mu(3).*x.^2 -1i*mu(2).*x + mu(1));
    
figure(2)
PhasePlot(z,F(z),'e'); % or e
h = gca;
set(h,'position',[0 0 1 1]);

print(['../../../general-presentations/images/animated-images/kernAnim',num2str(l)],'-dpng','-r50');
l
end