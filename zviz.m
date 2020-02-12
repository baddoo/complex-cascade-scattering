function zviz(varargin);
%ZVIZ  Visualize functions of a complex variable.
%   ZVIZ expr1 expr2 ... exprn or ZVIZ('expr1','expr2',...'exprn') plots
%   colored contour maps of the given expressions, which should involve
%   a complex variable, z.  Contour lines for the real and imaginary parts
%   are superimposed on a full color image using a color scheme developed
%   by John Richardson.  (http://home1.gte.net/jrsr/complex.html)
%   Red is real, blue is positive imaginary, green is negative imaginary,
%   black is small magnitude and white is large magnitude.  Branch cuts
%   appear as color discontinuities and coalescent contour lines.
%   Singularities appear as characteristic black and white patterns.
%
%   For example
%      zviz  z  2./z.^2  sqrt(z)  acosh(z)  atanh(z)  exp(-1./z.^2)
%   produces six figures with plots of the six expressions.
%
%   Regarding branch cuts, see pp. 79 and 86 of Abramowitz and Stegun
%   and Cleve's Corner, MATLAB News & Notes, Summer, 1998.
%
%                                   z = x + i*y
%   sqrt(z), log(z):                x <= 1       &  y == 0
%   asin(z), acos(z), atanh(z):     abs(x) >= 1  &  y == 0
%   atan(z), asinh(z):              x == 0       &  abs(y) >= 1
%   acsc(z), asec(z), acoth(z):     abs(x) <= 1  &  y == 0
%   acot(z), acsch(z):              x == 0       &  abs(y) <= 1
%   acosh(z):                       x <= 1       &  y == 0
%   asech(z):                  (x <= 0 | x >= 1) &  y == 0

warnings = warning;
warning('off')
fig = gcf;
for k = 1:nargin

   F = varargin{k};
   [x,y] = meshgrid(-2:1/32:2);
   z = x + i*y;
   w = feval(inline(F),z);

   figure(fig+k-1)
   image(z2rgb(w))
   title(F,'fontsize',12)
   axis square
   tic = 1:32:129;
   labels = ['-2';'-1';' 0';' 1';' 2'];
   set(gca,'xtick',tic,'xticklab',labels,'ytick',tic,'yticklab',labels)

   hold on
   c = -pi:pi/16:pi;
   contour(real(w),c,'k')
   contour(imag(w),c,'k')
   hold off

   drawnow
end
warning(warnings)


function C = z2rgb(z);
%Z2RGB Complex variable color image.
%   Z2RGB(Z) maps a complex matrix into a full color image.
%   The mapping is due to John Richardson
%   See http://home1.gte.net/jrsr/complex.html.
%   Example:
%      [x,y] = meshgrid(-2:1/16:2); z = x+i*y;
%      image(z2rgb(F(z)))

r = abs(z);
a = sqrt(1/6)*real(z);
b = sqrt(1/2)*imag(z);
d = 1./(1+r.^2);
R = 1/2 + sqrt(2/3)*real(z).*d;
G = 1/2 - d.*(a-b);
B = 1/2 - d.*(a+b);
d = 1/2 - r.*d;
d(r<1) = -d(r<1);
C(:,:,1) = R + d;
C(:,:,2) = G + d;
C(:,:,3) = B + d;
