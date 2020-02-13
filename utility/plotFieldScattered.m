function plotFieldScattered(h,ADData,AAData,plotData,type)

nP=12;

Beta = ADData.Beta; % PG factor
semiChordDim = ADData.chordDim/2;
% Put back into dimensional variables
X = plotData.X*semiChordDim;
Y = plotData.Y*semiChordDim./Beta; % Undo PG transformation.
sDim = ADData.spacDim(1);
dDim = ADData.spacDim(2);

Yper = funPer(Y.',sDim,nP);
Xper = funPer(X.',dDim,nP);

kx = AAData.kx;
ky = AAData.ky;
omega = AAData.omega;
M = ADData.M;

GM0 = -kx/Beta^2;
PM0 = -omega/Beta^2;

if strcmp(type,'acoustic')
    pangle = AAData.sigma;
elseif strcmp(type,'hvelocity')
    pangle = AAData.sigmao;
elseif strcmp(type,'vvelocity')
    pangle = AAData.sigmao;
elseif strcmp(type,'pressure')
    pangle = AAData.sigmao;
elseif strcmp(type,'source')
    pangle = AAData.sigma;
elseif strcmp(type,'periodic')
    pangle = 0;
end

h2 = funPerMult(h.',exp(1i*pangle),nP);

rot=exp(0i*atan(dDim/sDim));
newcoord = (Xper+1i*Yper)*rot;

h=pcolor(real(newcoord),imag(newcoord),real(h2));
set(h, 'EdgeColor', 'none');

hold on 

shading interp; colormap jet;

if isfield(plotData,'colLimits'); caxis(plotData.colLimits); end
if isfield(plotData,'axisLimits'); axis(plotData.axisLimits); end
if isfield(plotData,'colorbar'); colorbar; end

for l = -nP:nP
    loc = ([0,2*semiChordDim]+l*dDim)+1i*l*[sDim,sDim];
   plot(real(rot*loc),imag(rot*loc),'k','LineWidth',3) 
end

xlims = [-2*real(rot),3*real(rot)];
xlim(xlims);
axis equal
ax = gca;
ylims = ax.YLim;
drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'MaxHeadSize',1, varargin{:} );       
b = omega*ky*Beta; a = (-GM0+M^2*PM0);
l = rot.*(a+1i*b);
drawArrow(.99*xlims(1)+[0,sqrt(1)*cos(angle(l))],.9*ylims(1) + [0,sqrt(1)*sin(angle(l))],'linewidth',3,'color','k')

shading interp
hold off

xlabel('$x^\ast$','interpreter','LaTeX')
ylabel('$y^\ast$','interpreter','LaTeX')
%xlim([-1,2])
%real(rot)
%imag(rot)

%axis([-4*real(rot),5*real(rot),-4*imag(rot),4*imag(rot)])
%xlim([-2,2])

colorbar
%-ADData.M.*ADData.c0
%caxis(2e-1*[-1,1])

end