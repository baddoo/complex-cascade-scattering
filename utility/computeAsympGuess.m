% Calculates the asymptotic approximations for the roots of the
% Wiener--Hopf kernel function.

function asympGuess = computeAsympGuess(ADData,AAData,Modes)

% Extract relevant data from structures
mu = ADData.mu;
s = ADData.spac(1); d = ADData.spac(2);
del = sqrt(s^2+d^2);
chie = ADData.chie;
sigma = AAData.sigma;
n= 1:Modes.trunc;
expMC = exp(-1i*chie); expPC = exp(1i*chie);

% The Q suffix corresponds to the quadrant that the mode is in
switch ADData.case

    case 1
    a0q1  = 2i*pi*expMC/del;
    a1q1  = -expMC/del;
    a2q1  = (1i*(pi/2 + chie - sigma) - log(2*pi/mu(1)/del))*expMC./del;
    
    a0q2 = 2i*pi*expPC/del;
    a1q2 =  expPC/del;
    a2q2 = (1i*(pi/2 + chie + sigma) + log(2*pi/mu(1)/del))*expPC./del;
    
    a2a3 = (-1i*(pi/2 + chie + sigma) + log(2*pi/mu(1)/del))*expMC./del;
    a2q4 = (-1i*(pi/2 + chie - sigma) - log(2*pi/mu(1)/del))*expPC./del;
    
    case 2

    a0q1  = 2i*pi*expMC/del;
    a1q1  = 0;
    a2q1 = (-1i*sigma - log(1+1./(1i*mu(2))))*expMC./del;

    a0q2 = 2i*pi*expPC/del;
    a1q2 =  0;
    a2q2 = ( 1i*sigma + log(1-1./(1i*mu(2))))*expPC./del;
    
    a2a3 = (-1i*sigma + log(1-1./(1i*mu(2))))*expMC./del;
    a2q4 = ( 1i*sigma - log(1+1./(1i*mu(2))))*expPC./del;

    case 3
    a0q1  = 2i*pi*expMC/del;
    a1q1  = 0;
    a2q1  = -sigma*1i*expMC./del;

    a0q2 = 2i*pi*expPC/del;
    a1q2 =  0;
    a2q2 = sigma*1i*expPC./del;
    
    a2a3 = conj(a2q2);
    a2q4 = conj(a2q1);

end

asympQ1 = a0q1*n+a1q1*log(n)+a2q1;
asympQ2 = a0q2*n+a1q2*log(n)+a2q2;
asympQ3 = conj(a0q2)*n+conj(a1q2)*log(n)+a2a3;
asympQ4 = conj(a0q1)*n+conj(a1q1)*log(n)+a2q4;

asympGuess = [asympQ1,asympQ2,asympQ3,asympQ4];

end