function PPOnSphere(z,w,c,pres,tjmp,hght)
% phase plot of function on Riemann shphere
%
% Usage: PPOnSphere(z,w,c,pres,tjmp,hght)
%
% z - values on domain
% w - values of function at points z
% c - color scheme (optional)
% pres - resolution of phase (optional)
% tjmp -  jumps of phase (optional)
% hght -  cryptic scaling of radius to produce landscapes (optional)

% Part of the Complex Function Explorer (former Phase Plot) package 
% Version 1.1, February 1, 2014
% Copyright(c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set lght = 1 for more realistic representation with light
lght = 0;

if nargin<6, hght = 0; end

if nargin>=2 
  % convert scalars to constant function  
  if size(w)==1, w=w*ones(size(z)); end
    
  % change orientation of color scheme for correct zero-pole representation
  w = conj(w);
end

if nargin==6 
  if max(max(hght))<0, hght = -hght; end
  RGB = colscheme(w,c,tjmp,pres); 
  
elseif nargin==5 
  RGB = colscheme(w,c,tjmp,pres); 
  
elseif nargin==4
  RGB = colscheme(w,c,[],pres); 
  
elseif nargin==3
  RGB = colscheme(w,c,[],30); 
  
elseif nargin==2 
  RGB = colscheme(w,'m'); 
  
elseif nargin<=1 
  disp(' ')
  disp('The best way to make a phase plot on the entire sphere is')
  disp('the usage of two different charts (no singularities)')
  disp(' ')
  disp('First phase plot on lower hemishpere from transformed unit square:');
  disp(' ')
  disp('z = zdomain(-1-1i,1+1i,1000,1000);')
  disp('w = (z-1)./(z.^2+z+1);')
  disp('figure(1)')
  disp('hold on')
  disp('PPOnSphere(z,w);')

  z = zdomain(-1-1i,1+1i,1000,1000);
  w = (z-1)./(z.^2+z+1);
  
  figure(1)
  clf
  hold on
  PPOnSphere(z,w);
  axis(0.85*[-1,1,-1,1,-1,1])
  disp('paused, press key to continue');
  pause
  
  disp(' ')
  disp('Next phase plot on upper hemishpere from inverted unit square:');
  disp(' ')
  disp('z = 1./z;')
  disp('w = (z-1)./(z.^2+z+1);')
  disp('PPOnSphere(z,w);')
  z = 1./z;
  w = (z-1)./(z.^2+z+1);
  PPOnSphere(z,w);
  disp('Phase Plot on upper hemishpere added');
  disp(' ')
  axis(0.85*[-1,1,-1,1,-1,1])
  
  return
end

[p,q,r] = stereoP2S(z);

if length(hght)==1 % hght is scalar factor
  % stretch radius depending on modulus of function to get nice shapes
  % of analytic landscapes above sphere
  
  % spiky peaks and valleys
  HghtW = abs(w);
  
  % corresponds to (shifted) z-coordinate on Riemann sphere
  %HghtW = abs(w).^2;
  
  % flat peaks and valleys 
  %HghtW = abs(w).^3;
  
  % simple variant
  %sc = 1+hght*(HghtW-1)./(HghtW+1);
      
  % more sophisticated variant
  sc = (1/2)*(1+tan((pi/4)*(1+hght*(HghtW-1)./(HghtW+1))));
  
  h = surf(sc.*p,sc.*q,sc.*r,RGB);

elseif length(hght)>1 % hght is assumed to be precomputed scaling factor
  h = surf(hght.*p,hght.*q,hght.*r,RGB);
      
else
  h=surf(p,q,r,RGB);
end

set(h,'EdgeColor','none');

axis equal
axis off
view(40,20)

axis tight
setaxis=axis;
axis(0.85*setaxis);

if lght==1
    
lightangle(140,60)

set(gcf,'Renderer','zbuffer')
set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...
    'AmbientStrength',.7,'DiffuseStrength',.3,...
    'SpecularStrength',.5,'SpecularExponent',25,...
    'BackFaceLighting','unlit')

end

end

