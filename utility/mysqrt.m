function root = mysqrt(a,b)
%MYSQRT Calculates sqrt(a^2-b^2) with correct branch cut.
root1=sqrt(bsxfun(@minus,a,b).*exp(-1i*pi/2))/(exp(-1i*pi/4));
root2=sqrt(bsxfun(@plus,a,b).*exp(-1i*pi/2))/(exp(-1i*pi/4));
root=root1.*root2;
end

