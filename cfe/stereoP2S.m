function [p,q,r] = stereoP2S(z)
%  stereographic projection from plane to sphere
%
% Usage: [p,q,r] = stereoP2S(z)

% Part of the Complex Function Explorer (former Phase Plot) package 
% Version 1.1, February 1, 2014
% Copyright(c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = real(z);
y = imag(z);

w = x.^2 + y.^2 + 1;

p = 2*x./w;
q = 2*y./w;
r = (w-2)./w;
