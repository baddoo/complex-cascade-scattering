function RGB = brightenRGB(RGB,bright)
% modification of color scheme (internal use) 

% Usage: RGB = brightenRGB(RGB,bright)
%
% bright - scalar or field with values between 0 and 1 for brightening

% Part of the Complex Function Explorer (former Phase Plot) package 
% Version 1.1, February 1, 2014
% Copyright(c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(size(bright))==1
    bright = bright * ones(size(RGB(:,:,1)));
end
    
RGB(:,:,1) = (1-bright).*RGB(:,:,1) + bright.*ones(size(RGB(:,:,1))); 
   %+ (bright<0).*((1+bright).*RGB(:,:,2));
  

RGB(:,:,2) = (1-bright).*RGB(:,:,2) + bright.*ones(size(RGB(:,:,2)));
   %+ (bright<0).*((1+bright).*RGB(:,:,2));
   

RGB(:,:,3) = (1-bright).*RGB(:,:,3) + bright.*ones(size(RGB(:,:,3)));
    %+ (bright<0).*((1+bright).*RGB(:,:,3));

