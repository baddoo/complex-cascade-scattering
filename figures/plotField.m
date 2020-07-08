function plotField(fun,ADData,AAData,plotData,type)

nP=6;

Beta = ADData.Beta; % PG factor
semiChordDim = ADData.chordDim/2;

% Put back into dimensional variables
X = plotData.X*semiChordDim;
Y = plotData.Y*semiChordDim./Beta; % Undo PG transformation.
sDim = ADData.spacDim(1);
dDim = ADData.spacDim(2);

s = ADData.spac(1);
d = ADData.spac(2);

Yper = funPer(Y.',sDim,nP);
Xper = funPer(X.',dDim,nP);

kx = AAData.kx;
ky = AAData.ky;
omega = AAData.omega;
M = ADData.M;

GM0 = -real(kx)/Beta^2;
PM0 = -real(omega)/Beta^2;

funInc = exp(1i*(-GM0.*X/semiChordDim + omega*ky*Y/semiChordDim*Beta));

if strcmp(type,'potential')
    pangle = AAData.Sigma;
    funIncFin = funInc;
elseif strcmp(type,'hvelocity')
    pangle = AAData.Sigmao;
    funIncFin = (-1i*GM0+1i*M^2*PM0).*funInc.*exp(1i*M^2*PM0.*X/semiChordDim);
elseif strcmp(type,'vvelocity')
    pangle = AAData.Sigmao;
    funIncFin = (1i*omega*ky).*funInc.*exp(1i*M^2*PM0.*X/semiChordDim);
elseif strcmp(type,'pressure')
    pangle = AAData.Sigmao;
    funIncFin = 1i*(GM0 - PM0).*funInc.*exp(1i*M^2*PM0.*X/semiChordDim);
elseif strcmp(type,'source')
    pangle = AAData.Sigma;
elseif strcmp(type,'periodic')
    pangle = 0;
end

fun2 = funPerMult(fun(plotData.X+1i*plotData.Y).',exp(1i*pangle),nP);
hIncFinPer = funPerMult(funIncFin.',exp(1i*pangle),nP);

rot=exp(1i*atan(dDim/sDim));
newcoord = (Xper+1i*Yper)*rot;

h=pcolor(real(newcoord),imag(newcoord),real(fun2 + hIncFinPer));
set(h, 'EdgeColor', 'none');

hold on 

%pcolor(real(newcoord),imag(newcoord),0*real(hIncFinPer));
shading interp; colormap jet;

if isfield(plotData,'colLimits'); caxis(plotData.colLimits); end
if isfield(plotData,'axisLimits'); axis(plotData.axisLimits); end
if isfield(plotData,'colorbar'); colorbar; end


for l = -nP:nP
    loc = ([0,2*semiChordDim]+l*dDim)+1i*l*[sDim,sDim];
   plot(real(rot*loc),imag(rot*loc),'k','LineWidth',10) 
end

axis equal

xlims = 2*semiChordDim*[-2*real(rot),3*real(rot)];
xlim(xlims);

ax = gca;
ylims = ax.YLim;
drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'MaxHeadSize',1, varargin{:} );       
b = omega*ky*Beta; a = (-GM0+M^2*PM0);
l = rot.*(a+1i*b);
%drawArrow(.95*xlims(1)+[0,sqrt(.5)*cos(angle(l))],.95*ylims(2) + [0,sqrt(.5)*sin(angle(l))],'linewidth',3,'color','k')
start = [.95*xlims(1),.95*ylims(2)];
stop = [.95*xlims(1)+sqrt(.5)*cos(angle(l)),(.95*ylims(2) + sqrt(.5)*sin(angle(l)))];
shading interp
arrow(start,stop,30,'LineWidth',8);

hold off

xlabel('$x/2$','interpreter','LaTeX')
ylabel('$y/2$','interpreter','LaTeX')

colorbar

end