function  z = zplanePP(n)
% discrete field of points covering the entire complex plane
%
% Usage: z = zplanePP(n);

% Part of the Complex Function Explorer (former Phase Plot) package 
% Version 1.1, February 1, 2014
% Copyright(c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1, n=1000; end;

[p,q,r] = sphere(n);

z = stereoS2P(p,q,r);

