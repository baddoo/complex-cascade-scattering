function PP = pltphase(z,f,cs,t,pres,height)
% phase plot of complex function f(z), internal use 

% Usage: PP = pltphase(z,f,cs,t,pres,height)

% Part of the Complex Function Explorer (former Phase Plot) package 
% Version 1.1, February 1, 2014
% Copyright(c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6, height=0; end

if nargin<5 || isempty(pres), pres=20; end

if nargin<4 || isempty(t), t=exp(1i*linspace(-pi,pi,pres+1)); t=t(1:pres); end

if nargin<3 || isempty(cs), cs = 'p'; end

if nargin<2, z=zdomain; f=sin(5*z); end

% convert scalar to a constant function  
if length(f)==1, f=f*ones(size(z)); end

if nargin<6
  % if height is not given, choose shape of surface on which the plot lives

  %rep='a';  standard analytic landscape 
  %rep='l';  logarithmic analytic landscape 
  %rep='c';  compressed analytic landscape 
  
  rep='0';  height=zeros(size(z)); hmin=0; hmax=1; %#ok<NASGU>
     
  if length(cs)>=2
  
    % modify if color scheme requests special surface
    incsa = strfind(cs(2:end),'a');
    if ~isempty(incsa)>0, rep=cs(incsa(1)+1); cs(incsa+1) = ''; end

    incsl = strfind(cs(2:end),'l');
    if ~isempty(incsl)>0, rep=cs(incsl(1)+1); cs(incsl+1) = ''; end

    incsc = strfind(cs(2:end),'c');
    if ~isempty(incsc)>0, rep=cs(incsc(1)+1); cs(incsc+1) = ''; end

    if rep=='a'
      % analytic landscape, cutoff at hmax
      hmax=5;
      height=abs(f); 
      height(height>hmax)=hmax;
      % normalize to [0,1]
      hmin = min(min(height)); hmax = max(max(height));
      if hmax>hmin
        height = (1/(hmax-hmin))*(height-hmin);
      end
      
    elseif rep=='l'
      % logarithmic analytic landscape
      hmin =-3; hmax=3;
      absf = abs(f);
      height = log(absf);
      small = (height<hmin);
      large = (height>hmax);
      height(small) = hmin;
      height(large) = hmax;
      % normalize height to [0,1]
      hmin = min(min(height)); hmax = max(max(height));
      if hmax>hmin
        height = (1/(hmax-hmin))*(height-hmin);
      end
      
    elseif rep=='c' 
      % compressed analytic landscape, zeros and poles at levels 0 and 1
      height=(2/pi)*atan(abs(f));
    
    end
   
  end
  
end

% flip arrays for correct orientation of surface normal 
z=flipud(z);
f=flipud(f);
height=flipud(height);

hmin = min(min(height));
hmax = max(max(height));

% scale height for nice shape of surface
xmin=min(min(real(z)));
xmax=max(max(real(z)));
ymin=min(min(imag(z)));
ymax=max(max(imag(z)));

scale=0.4*(xmax-xmin+ymax-ymin);
height=scale*height;

% get RGB colors according to chosen color scheme
RGB  = colscheme(f,cs,t,pres); 

% plot surface
PP = surf(real(z),imag(z),height,RGB);
set(PP,'EdgeColor','none');

% set appropriate axis
hmin=scale*hmin;
hmax=scale*hmax;

axis equal
axis tight
ax = axis;
if hmax>hmin, axis([ax,hmin,hmax]); end
axis off


