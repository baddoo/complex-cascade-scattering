% The logarithmic derivative of numerator of K.
% The domain must be partitioned, and different realisations of K applied
% to avoid rounding errors.

function KLD = KNumlogD(gamma,ADData,AAData)

% Make a zero matrix of the correct size
KLD = zeros(size(gamma));

% Extract data from structures
omega = AAData.omega;
w = AAData.w;
he = ADData.spac(1);
d = ADData.spac(2);
Sigma = AAData.Sigma;
mu = ADData.mu; 

zeta = @(x) mysqrt(omega*w,x);
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

zg1 = zg([quad1(:);quad3(:)]);
zg2 = zg(quad2);
zg4 = zg(quad4);

g1 = gamma([quad1(:);quad3(:)]);
g2 = gamma(quad2);
g4 = gamma(quad4);

% Now we define all the functions in their respective regions
    
% Define polynomial
p1 = mu(1) - 1i*g1*mu(2) - g1.^2*mu(3);
p2 = mu(1) - 1i*g2*mu(2) - g2.^2*mu(3);
p4 = mu(1) - 1i*g4*mu(2) - g4.^2*mu(3);

% Define deriative of polynomial
Dp1 = -1i*mu(2) - 2*g1*mu(3);
Dp2 = -1i*mu(2) - 2*g2*mu(3);
Dp4 = -1i*mu(2) - 2*g4*mu(3);

expPHPD1 = exp(1i*(he*zg1 + d*g1));
expPHPD2 = exp(1i*(he*zg2 + d*g2));
expPHPD4 = exp(1i*(he*zg4 + d*g4));

expPHND1 = exp(1i*(he*zg1 - d*g1));
expPHND2 = exp(1i*(he*zg2 - d*g2));
expPHND4 = exp(1i*(he*zg4 - d*g4));

expPH1 = exp(1i*he*zg1);
%expPH2 = exp(1i*he*zg2);
%expPH4 = exp(1i*he*zg4);

%expPD1 = exp(1i*d*g1);
expPD2 = exp(1i*d*g2);
expPD4 = exp(1i*d*g4);

expPS = exp(1i*Sigma);

% Define regularisations of denominator
den1 = expPH1.^2+1 - expPHPD1*expPS - expPHND1/expPS;
den2 = expPHPD2 + 1./expPHND2 - expPD2.^2*expPS - 1./expPS;
den4 = expPHND4 + 1./expPHPD4 - expPS - expPD4.^-2/expPS;

% Define regularisations of first part of numerator
num1 = zg1/1i.*(expPH1.^2-1);
num2 = zg2/1i.*(expPHPD2-1./expPHND2);
num4 = zg4/1i.*(expPHND4-1./expPHPD4);

% Define regularisations of derivative of first part of denominator
dNum1R = -g1./zg1/1i.*(expPH1.^2-1) - he*g1.*(expPH1.^2+1);
dNum2  = -g2./zg2/1i.*(expPHPD2 - 1./expPHND2) - he*g2.*(expPHPD2 + 1./expPHND2);
dNum4  = -g4./zg4/1i.*(expPHND4 - 1./expPHPD4) - he*g4.*(expPHND4 + 1./expPHPD4);

dDenR1 = he*g1./zg1/1i.*(expPH1.^2-1) + d/1i*(expPHPD1*expPS - expPHND1/expPS);
dDenR2 = he*g2./zg2/1i.*(expPHPD2 - 1./expPHND2) + d/1i*(expPD2.^2*expPS- 1/expPS);
dDenR4 = he*g4./zg4/1i.*(expPHND4 - 1./expPHPD4) + d/1i*(expPS - expPD4.^-2/expPS);

fullNum1 = num1 + p1.*den1;
fullNum2 = num2 + p2.*den2;
fullNum4 = num4 + p4.*den4;

fullDNum1 = dNum1R + p1.*dDenR1 + Dp1.*den1;
fullDNum2 = dNum2 +  p2.*dDenR2 + Dp2.*den2;
fullDNum4 = dNum4 +  p4.*dDenR4 + Dp4.*den4;

% Fill in quadrants
KLD([quad1(:);quad3(:)]) = fullDNum1./fullNum1;
KLD(quad2) = fullDNum2./fullNum2;
KLD(quad4) = fullDNum4./fullNum4;

end

end

% Elapsed time is 1.762890 seconds.
% Elapsed time is 1.095508 seconds.
% Elapsed time is 1.189526 seconds.


