function [asympRoots,asympGuess,boxHeight] = findAsympRoots(ADData,AAData,Modes)

mu = ADData.mu;
he = ADData.spac(1);
d = ADData.spac(2);
sigma5 = AAData.sigma5;
w = AAData.w;
omega5 = AAData.omega5;

sols1ab = [];
asympGuess = computeAsympGuess(ADData,AAData,Modes);

for nL = 1:numel(asympGuess)
sols1ab = [sols1ab,findRoots3(he,d,sigma5,mu,omega5,w,asympGuess(nL))];
end

asympRoots = sols1ab;

position = asympGuess(real(asympGuess)<0);
boxHeight = 6*abs(imag(position(1)));

end