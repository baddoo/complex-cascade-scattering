function Kminus = Kminus(gammav,Kargs)
%% Extract data from structs
TP=Kargs.TP; LP=Kargs.LP;

%% Define numerator
num=1-bsxfun(@rdivide,gammav,TP);
numprod=prod(num,3);

%% Define denominator
den=1-bsxfun(@rdivide,gammav,LP);
denprod=prod(den,3);

%% Define constant and exponential term
Kminus=numprod./denprod;

end