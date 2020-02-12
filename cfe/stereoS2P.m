function zz = stereoS2P(x,y,z)
%  stereographic projection from sphere to plane
%
% Usage: zz = stereoS2P(x,y,z)

% Part of the Complex Function Explorer (former Phase Plot) package 
% Version 1.1, February 1, 2014
% Copyright(c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zz = (x+1i*y)./(1-z);
