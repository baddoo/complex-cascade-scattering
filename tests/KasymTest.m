addpath(genpath('../'))

imageFolder = '../../images/';

chordDim = 1;
stagAngDim = 50;

vaneDistDim = chordDim*.2;

dDim = vaneDistDim*sin(stagAngDim*pi/180);
sDim = vaneDistDim*cos(stagAngDim*pi/180);

%l
%%
ADData=struct('spacDim',     [sDim,dDim],... %Dimensional blade spacing
              'chordDim',    chordDim,...           %Blade length
              'c0',          340,...         %Speed of sound
              'M',           0.4,...           %Mach number 
              'Wdim',        0,...           %Spanwise background velocity
              'case',        1, ...              %Case
              'C',           .5);                 %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omegaDim',    [],...                 %Time Frequency
               'omega',       1,...                 %Time Frequency
               'kxDim',       [],...                   %Tangential frequency
               'kx',          1,...                   %Tangential frequency
               'ky',          [],...
               'kz',          [],...
               'kyDim',       [],...                       %Normal frequency
               'kzDim',       0,...                       %Spanwise frequency
               'sigmao',      3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Amp',         [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',20000,...                     %Truncation of kernel modes
             'dmodes',5,...                      %Number of duct mode
             'amodes',5);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);

%%
[asympGuess,coefs] = computeAsympGuess(newADData,newAAData,Modes);
largeGammaP =  1i*logspace(0,4,1e2);
KlargeP = Kplus(largeGammaP,out.Kargs);
largeGammaN = -1i*logspace(0,4,1e2);
KlargeN = Kminus(largeGammaN,out.Kargs);
mu = newADData.mu;
cse = newADData.case; chie = newADData.chie;

switch cse
    case 1
        expPI = .5; + chie/pi;
        expNI = .5 ;+ chie/pi;
    case 2
        % Need to choose branch here correctly.
        expPI =  1/pi*acot(mu(2))+1;%+(1+-sign(-angle(1i-mu(2))))/2;
        expNI = -1/pi*acot(mu(2));%+(1-sign(-angle(1i-mu(2))))/2;
    case 3
        expPI = 1;
        expNI = 1;
end

%expPI = +coefs.a2q1./coefs.a0q1 + coefs.a2q2./coefs.a0q2;
%expNI = +coefs.a2q3./coefs.a0q3 + coefs.a2q4./coefs.a0q4;

% exponentP =mod(real(expPI)-2*eps,1)+1i*imag(expPI);
% exponentN = mod(real(expNI)-2*eps,1)+1i*imag(expNI);
% exponentN = 2+expNI;

KasympP = K(0,newADData,newAAData)*exp((expPI)*log(largeGammaP)); % Asymptotic approximation up to a multiplativative constant
KasympN = exp((expNI)*log(largeGammaN)); % Asymptotic approximation up to a multiplativative constant

% Final result
figure(1)
subplot(2,1,1);
loglog(abs(largeGammaP),abs(KlargeP),'b','LineWidth',3);
hold on
loglog(abs(largeGammaP),abs(KasympP),'b--','LineWidth',5);
loglog(abs(largeGammaN),abs(KlargeN),'r','LineWidth',3);
loglog(abs(largeGammaN),abs(KasympN),'r--','LineWidth',3);
axis square
hold off

subplot(2,1,2);
loglog(abs(largeGammaP),abs(KlargeP./KasympP),'b','LineWidth',3);
hold on
loglog(abs(largeGammaN),abs(KlargeN./KasympN),'r','LineWidth',3);
axis square
hold off


%% Verification of proposition 2
largeGammaM = -(1+1i)*logspace(00,7,1e3);
TPR = out.TP(real(out.TP)>=0); TPR = TPR(1:3000);
TPL = out.TP(real(out.TP)<0); TPL = TPL(1:3000);
mR3 = permute(1:numel(TPR),[1,3,2]);
mL3 = permute(1:numel(TPL),[1,3,2]);
asymR = coefs.a0q1 * (mR3-1) + 0*coefs.a1q1 * log(mR3) + coefs.a2q1;
asymL = coefs.a0q2 * (mL3-1) + 0*coefs.a1q2 * log(mL3) + coefs.a2q2;
%asymL = coefs.a0q2 * mL3 + coefs.a2q2;
%asymR = coefs.a0q1 * mR3 + coefs.a2q1;

P2 = prod(1./(TPR-largeGammaM).*(asymR-largeGammaM)./(TPL-largeGammaM).*(asymL-largeGammaM),3);
%P2 = prod
figure(2)
%semilogy(abs(permute(TPR-asymR,[1,3,2])))
loglog(abs(largeGammaM),abs(P2-1),'b','LineWidth',3);


