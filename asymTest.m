n3 = permute(1:5e3,[1,3,2]);
a = 1.5; b = .85; c = 0; d =2;
an3 = a*n3 + b + c*log(n3) + d./n3;
nz = 100;
z = linspace(0,100,nz)*exp(-1i*pi/4);
exact1 = prod((1-bsxfun(@rdivide,z,an3)).*exp(bsxfun(@rdivide,z,an3)),3);
%exact2 = prod((1-z./(a*n3+b)),3);
approx = exp((-z./a+b./a+.5).*log(-a./z)-z./a*(1-eulergamma)+z.*sum(1./an3-1./(a*n3),3));
%approx = exp(-(-z./a+b./a+.5).*log(-z./a)-z./a*(1-eulergamma)+z.*sum(-(b+d./n3)./(a*n3.*(a*n3+b+d./n3)),3));

%myProd = prod((an3-z)./(a*n3+b+c*log(n3)-z),3);

%semilogy(abs(z),abs(exact1))
plot(abs(exact1./approx))
hold on
%semilogy(abs(approx))
hold off
