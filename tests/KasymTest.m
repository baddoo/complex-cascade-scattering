imageFolder = '../../images/';

chordDim = .6;
stagAngDim = 40;

vaneDistDim = chordDim/cos(stagAngDim*pi/180);

dDim = vaneDistDim*sin(stagAngDim*pi/180);
sDim = vaneDistDim*cos(stagAngDim*pi/180);

%l
%%
ADData=struct('spacDim',     [sDim,dDim],... %Dimensional blade spacing
              'chordDim',    chordDim,...           %Blade length
              'c0',          340,...         %Speed of sound
              'M',           0.5,...           %Mach number 
              'Wdim',        0,...           %Spanwise background velocity
              'case',        2, ...              %Case
              'C',           -2i);                 %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omegaDim',    [],...                 %Time Frequency
               'omega',       30,...                 %Time Frequency
               'kxDim',       [],...                   %Tangential frequency
               'kx',          5,...                   %Tangential frequency
               'kyDim',       [],...                       %Normal frequency
               'kzDim',       0,...                       %Spanwise frequency
               'sigmao',      3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Amp',         [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',2000,...                     %Truncation of kernel modes
             'dmodes',5,...                      %Number of duct mode
             'amodes',5);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);

%%
[asympGuess,coefs] = computeAsympGuess(newADData,newAAData,Modes);
largeGammaP = 1i*logspace(0,3,1e2);
KlargeP = Kplus(largeGammaP,out.Kargs);
largeGammaN = -1i*logspace(0,3,1e2);
KlargeN = Kminus(largeGammaN,out.Kargs);
mu = newADData.mu;
cse = newADData.case; chie = newADData.chie;
switch cse
    case 1
        exponentPI = .5 + chie/pi;
        exponentNI = .5 + chie/pi;
    case 2
        exponentPI = 1/pi*acot(mu(2))+1;
        exponentNI = 1/pi*acot(mu(2));
    case 3
        exponentPI = 1;
        exponentNI = 1;
end

%exponenta = coefs.a2q1./coefs.a0q1 + coefs.a2q2./coefs.a0q2;
%exponentb = coefs.a2q3./coefs.a0q3 + coefs.a2q4./coefs.a0q4;

exponentP =mod(real(exponentPI)-2*eps,1)+1i*imag(exponentPI);
exponentN = mod(real(exponentNI)-2*eps,1)+1i*imag(exponentNI);
exponentN = 2+exponentNI;

KasympP = exp((exponentPI)*log(largeGammaP)); % Asymptotic approximation up to a multiplativative constant
KasympN = exp((exponentNI)*log(largeGammaN)); % Asymptotic approximation up to a multiplativative constant

% Final result
figure(1)
subplot(2,1,1);
loglog(abs(largeGammaP),abs(KlargeP),'b','LineWidth',3);
hold on
loglog(abs(largeGammaP),abs(KasympP),'b--','LineWidth',5);
loglog(abs(largeGammaN),abs(KlargeN),'r','LineWidth',3);
loglog(abs(largeGammaN),abs(KasympN),'r--','LineWidth',3);
hold off

subplot(2,1,2);
loglog(abs(largeGammaP),abs(KlargeP./KasympP),'b','LineWidth',3);
hold on
loglog(abs(largeGammaN),abs(KlargeN./KasympN),'r','LineWidth',3);
hold off

return
%% Verification of proposition 2
largeGammaM = -1i*linspace(00,1e6,1e3);
TPR = out.TP(real(out.TP)>0);
TPL = out.TP(real(out.TP)<0);
mR3 = permute(1:numel(TPR),[1,3,2]);
mL3 = permute(1:numel(TPL),[1,3,2]);
asymR = coefs.a0q1*(mR3-2) + coefs.a1q1*log(mR3) + coefs.a2q1;
asymL = coefs.a0q2*(mL3+2) + coefs.a1q2*log(mL3) + coefs.a2q2;
P2 = prod((TPR-largeGammaM)./(asymR-largeGammaM),3)...
   .*prod((TPL-largeGammaM)./(asymL-largeGammaM),3);
%P2 = prod
figure(2)
%semilogy(abs(permute(TPR-asymR,[1,3,2])))
loglog(abs(largeGammaM),abs(P2),'b','LineWidth',3);


