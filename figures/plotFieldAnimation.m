function plotFieldTime(h,omeg,ADData,AAData,plotData,type)

nP=12;

Beta = ADData.Beta; % PG factor

X = plotData.X;
Y = plotData.Y./Beta; % Undo PG transformation.
dPhys = ADData.spac(2);
sPhys = ADData.spac(1)/Beta;

Yper = funPer(Y.',sPhys,nP);
Xper = funPer(X.',dPhys,nP);

delta = ADData.delta;
kx = AAData.k;
kn = AAData.kn;
omega = AAData.omega;
M = ADData.M;

GM0 = -delta.*kx;
PM0 = -delta.*omega;

hInc = exp(1i*(-GM0.*X+omega*kn*Y.*Beta));

if strcmp(type,'acoustic')
    pangle = AAData.sigma5;
    hIncFin = hInc;
elseif strcmp(type,'hvelocity')
    pangle = AAData.sigmao5;
    hIncFin = (-1i*GM0+1i*M^2*PM0).*hInc.*exp(1i*M^2*PM0.*X);
elseif strcmp(type,'vvelocity')
    pangle = AAData.sigmao5;
    hIncFin = (1i*omega*kn).*hInc.*exp(1i*M^2*PM0.*X);
elseif strcmp(type,'pressure')
    pangle = AAData.sigmao5;
    hIncFin = (-1i*GM0+1i*PM0).*hInc.*exp(1i*M^2*PM0.*X);
elseif strcmp(type,'source')
    pangle = AAData.sigma;
elseif strcmp(type,'periodic')
    pangle = 0;
end

h2 = funPerMult(h.',exp(1i*pangle),nP);
hIncFinPer = funPerMult(hIncFin.',exp(1i*pangle),nP);

rot=exp(1i*atan(dPhys/sPhys));
newcoord = (Xper+1i*Yper)*rot;

h=pcolor(real(newcoord),imag(newcoord),real(exp(-1i*omeg)*(h2 + hIncFinPer)));
set(h, 'EdgeColor', 'none');

hold on 

shading interp; colormap jet;

if isfield(plotData,'colLimits'); caxis(plotData.colLimits); end
if isfield(plotData,'axisLimits'); axis(plotData.axisLimits); end
if isfield(plotData,'colorbar'); colorbar; end


for l = -nP:nP
    loc = ([0,1]+l*dPhys)+1i*l*[sPhys,sPhys];
   plot(real(rot*loc),imag(rot*loc),'k','LineWidth',3) 
end

xlims = [-2*real(rot),3*real(rot)];
xlim(xlims);
axis equal
ax = gca;
ylims = ax.YLim;
drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'MaxHeadSize',1, varargin{:} );       
b = omega*kn*Beta; a = (-GM0+M^2*PM0);
l = rot.*(a+1i*b);
drawArrow(-1.3+.99*xlims(1)+[0,sqrt(1)*cos(angle(l))],.9*ylims(1) + [0,sqrt(1)*sin(angle(l))],'linewidth',3,'color','k')

hold off
axis off
set(gca,'position',[0 0 1 1])


end