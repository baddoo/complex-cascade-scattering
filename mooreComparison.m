
beta0 = 0;
beta1 = 1;
sig = 1;
C = @(sigVar) besselk(1,1i*sigVar)/(besselk(0,1i*sigVar) + besselk(1,1i*sigVar));

a = zeros(1,3);
a(1) = -2i*pi*C(sig)*beta0 + 2i*pi*(1-C(sig))*beta1 - 2*C(sig)*beta1;
a(2) = 2*pi^2*beta0 - 4i*pi*beta1;
a(3) = pi^2*beta1;

pMoore = @(xVar) a(1)*sqrt((1-xVar)./(1+xVar)) ...
          + a(2)*2*sqrt(1-xVar.^2) ...
          + a(3)*4*xVar.*sqrt(1-xVar.^2);
      
%% Plots

xPlot = sin(linspace(-1,1)*pi/2);
plot(xPlot,real(pMoore(xPlot)))