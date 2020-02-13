function [ root ] = mysqrt(a,b)
%UNTITLED Summary of this function goes here
% This function calculates the square root taking the correct branch cut. 
% The typical form will be a^2-b^2
% This assumes that the branch point has positive real part
root1=sqrt(bsxfun(@minus,a,b).*exp(-1i*pi/2))/(exp(-1i*pi/4));
root2=sqrt(bsxfun(@plus,a,b).*exp(-1i*pi/2))/(exp(-1i*pi/4));
root=root1.*root2;
end

