function w = MyFunction(z)
% example of complex function written as script for use with CFEGUI
%
% this is an atomic singular inner function with atoms sitting at the
% fifth root of unity

% Part of the Complex Function Explorer (former Phase Plot) package 
% Version 1.1, February 1, 2014
% Copyright(c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

omega = exp(2i*pi/5);

w = ones(size(z));

for k=0:4
  w = w.*exp((z+omega.^k)./(z-omega.^k));
end

end
