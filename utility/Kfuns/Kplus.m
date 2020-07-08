function Kplus = Kplus(gammav,Kargs)
% Factorisation of Wiener--Hopf kernel into the upper half plane.

%% Extract data from structs
TM=Kargs.TM; LM=Kargs.LM;

%% Define numerator
num=1-bsxfun(@rdivide,gammav,TM);
numprod=prod(num,3);
%% Define denominator
den = 1-bsxfun(@rdivide,gammav,LM);
denprod=prod(den,3);

%% Define constant and exponential term
const = K(0,Kargs.ADData,Kargs.AAData);
Kplus=const.*numprod./denprod;

end