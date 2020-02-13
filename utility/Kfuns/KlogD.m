% The logarithmic derivative of K.
% The domain must be partitioned, and different realisations of K applied
% to avoid rounding errors.

function KLD = KlogD(gamma,ADData,AAData)

% Make a zero matrix of the correct size
KLD = zeros(size(gamma));

% Extract data from structures
omega5 = AAData.omega5;
w = AAData.w;
he = ADData.spac(1);
d = ADData.spac(2);
sigma = AAData.sigma5;
mu = ADData.mu; 

zeta = @(x) mysqrt(omega5*w,x);
zg = zeta(gamma);

if mu == zeros(size(mu))
% In the special case where mu vanishes, we may write the logarithmic
% derivative in the simple, well conditioned form: 
KLD = -gamma.*(1+zg*he.*cot(he*zg))./zg.^2;   
    
else

% Otherwise, there are four regions to consider where different functions dominate and 
% there are difference scaling behaviours. Therefore we have to divide gamma 
% into quadrants. We label these 1 to 4 going anticlockwise.

quad1 = find(-atan(he/d)<=angle(gamma) & angle(gamma) <=atan(he/d));
quad2 = find(pi-atan(he/d)>=angle(gamma) & angle(gamma)>=atan(he/d));
quad3 = find(-atan(he/d)<=angle(-gamma) & angle(-gamma)<=+atan(he/d));
quad4 = find(-atan(he/d)>=angle(gamma) & angle(gamma)>=-pi+atan(he/d));

% Now we define all the functions in their respective regions
    
% Define polynomial
p = @(x) mu(1) + x*mu(2) + x.^2*mu(3);

% Define deriative of polynomial
Dp= @(x) mu(2) + 2*x*mu(3);

% Define regularisations of denominator
den1 = @(x) exp(2i*he*zeta(x))+1 - exp(1i*(he*zeta(x)+ d*x + sigma))-exp(-1i*(-he*zeta(x)+ d*x + sigma));
den2 = @(x) exp(1i*(he*zeta(x) + d*x))+exp(1i*( d*x-he*zeta(x))) - exp(1i*(2*d*x + sigma))-exp(-1i*(0*d*x + sigma));
den4 = @(x) exp(1i*(he*zeta(x) - d*x))+exp(1i*(-d*x-he*zeta(x))) - exp(1i*(0*d*x + sigma))-exp(-1i*(2*d*x + sigma));

% Define regularisations of first part of numerator
num1 = @(x) zeta(x)/1i.*(exp(2i*he*zeta(x))-1);
num2 = @(x) zeta(x)/1i.*(exp(1i*(he*zeta(x)+d*x))-exp(1i*( d*x-he*zeta(x))));
num4 = @(x) zeta(x)/1i.*(exp(1i*(he*zeta(x)-d*x))-exp(1i*(-d*x-he*zeta(x))));

% Define regularisations of derivative of first part of denominator
dNum1R = @(x) -x./zeta(x)/1i.*(exp(2i*he*zeta(x))-1) - he*x.*(exp(2i*he*zeta(x))+1);
dNum2  = @(x) -x./zeta(x)/1i.*(exp(1i*(he*zeta(x)+d*x))-exp(1i*( d*x-he*zeta(x)))) - he*x.*(exp(1i*(he*zeta(x)+d*x))+exp(1i*( d*x-he*zeta(x))));
dNum4  = @(x) -x./zeta(x)/1i.*(exp(1i*(he*zeta(x)-d*x))-exp(1i*(-d*x-he*zeta(x)))) - he*x.*(exp(1i*(he*zeta(x)-d*x))+exp(1i*(-d*x-he*zeta(x))));

dDenR1 = @(x) he*x./zeta(x)/1i.*(exp(2i*he*zeta(x))-1) + d/1i*(exp(1i*(he*zeta(x)+d*x + sigma))-exp(-1i*(-he*zeta(x)+d*x + sigma)));
dDenR2 = @(x) he*x./zeta(x)/1i.*(exp(1i*(he*zeta(x)+d*x))-exp(1i*( d*x-he*zeta(x)))) + d/1i*(exp(1i*(2*d*x + sigma))-exp(-1i*(        sigma)));
dDenR4 = @(x) he*x./zeta(x)/1i.*(exp(1i*(he*zeta(x)-d*x))-exp(1i*(-d*x-he*zeta(x)))) + d/1i*(exp(1i*(0*d*x + sigma))-exp(-1i*(2*d*x+ sigma)));

fullNum1 = @(x) num1(x) + p(x).*den1(x);
fullNum2 = @(x) num2(x) + p(x).*den2(x);
fullNum4 = @(x) num4(x) + p(x).*den4(x);

fullDNum1 = @(x) dNum1R(x) + p(x).*dDenR1(x) + Dp(x).*den1(x);
fullDNum2 = @(x) dNum2(x) + p(x).*dDenR2(x) + Dp(x).*den2(x);
fullDNum4 = @(x) dNum4(x) + p(x).*dDenR4(x) + Dp(x).*den4(x);

% Fill in quadrants
KLD(quad1) = fullDNum1(gamma(quad1))./fullNum1(gamma(quad1));
KLD(quad2) = fullDNum2(gamma(quad2))./fullNum2(gamma(quad2));
KLD(quad3) = fullDNum1(gamma(quad3))./fullNum1(gamma(quad3));
KLD(quad4) = fullDNum4(gamma(quad4))./fullNum4(gamma(quad4));

end

end
