%addpath(genpath('../'));
imageFolder = '../../../LaTeX/porous-jfm-r2s/images/';
addpath('../matlab2tikz/src')
chordDim = 1;
vaneDistDim = 0.6;
stagAngDim = 40;

nFreq = 100;
C = -[0.0001,0.01,.1,1];
nC = numel(C);
omegaLoop = linspace(0,30,nFreq+1); omegaLoop(1)=[];
lift = zeros(nC,nFreq);
myAModes = 35;
SPOU = zeros(nC,nFreq,2*myAModes+1);
SPOD = zeros(nC,nFreq,2*myAModes+1);

for m = 1:nC
for l = 1:nFreq
disp(['Loop number ' num2str(m) 'x' num2str(l)])
%%
ADData=struct('spacDim',    [vaneDistDim*cos(stagAngDim*pi/180),vaneDistDim*sin(stagAngDim*pi/180)],...%Blade Spacing
              'chordDim',   chordDim,...                          %Blade length
              'c0',         340,...                        %Speed of sound
              'M',          .3,...                         %Mach number 
              'Wdim',       0,...
              'case',       2,...
              'C',          C(m));                           
              
%% Aeroacoustic Data
AAData=struct( 'omegaDim', [],...
               'omega',    omegaLoop(l)+1e-6i,...                      %Tangential Frequency
               'kx',       'gust',...
               'kxDim',    [],...
               'ky',       [],...                       %Normal frequency
               'kyDim',    [],...               
               'kz',       0,...                       %Spanwise frequency
               'kzDim',    [],...                       %Spanwise frequency
               'Sigmao',   3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Sigma',    [],...
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',500,...                     %Truncation of kernel modes
             'dmodes',35,...                      %Number of duct mode
             'amodes',myAModes);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);    
data=computeCoefficients(newADData,newAAData,out);

%%
lift(m,l) = 1i*omegaLoop(l).*D(-omegaLoop(l).*newADData.M^2./newADData.Beta.^2,data)./newAAData.Amp(2);
LPa = out.LPa; LMa = out.LMa; ZPa = out.ZPa; ZMa = out.ZMa; SQRTa = out.SQRTa;

data.comb = [1 0 1 0];
DLP = permute(D(permute(LPa,[3,2,1]),data),[3,2,1]);
data.comb = [0 1 0 1];
DLM = permute(D(permute(LMa,[3,2,1]),data),[3,2,1]);

Beta = newADData.Beta;
s = sqrt(newADData.spac(1).^2+newADData.spac(2)^2);
se = newADData.spac(3);

spoU = pi^2*Beta^2*omegaLoop(l).*abs(ZPa.*DLP).^2.*real(1./SQRTa)./se./abs(newAAData.Amp(2).^2);
spoD = pi^2*Beta^2*omegaLoop(l).*abs(ZMa.*DLM).^2.*real(1./SQRTa)./se./abs(newAAData.Amp(2).^2);

SPOU(m,l,:) = spoU;
SPOD(m,l,:) = spoD;

end
end

%% Plot lift
figure(1)

c = hot(ceil(1.5*nC));
plot(omegaLoop,abs(lift(1,:)),'k','LineWidth',3);
    legendInfo{1} = ['$C_{II} = ' num2str(0) '$']; 
hold on
for j = 1:(nC-1)
plot(omegaLoop,abs(lift(j+1,:)),'LineWidth',2,'Color',c(j,:));
    legendInfo{j+1} = ['$C_{II} = ' num2str(C(j+1)) '$']; 
end

legend(legendInfo,'Interpreter','latex','Location','northeast')
hold off
ylim([0,1.5])
xlim([0,20])

xlabel('$\omega$');%,'Interpreter','Latex')
ylabel('$\left|C_L\right|$');%,'Interpreter','LaTeX')

hold off
cleanfigure
matlab2tikz([imageFolder,'ul-glegg','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

%% Plot first sound power output downstream mode
firstDSMode = SPOD(:,:,ceil(end/2));
secondDSMode = SPOD(:,:,ceil(end/2)+1);

figure(2)
plot(omegaLoop,real(firstDSMode(1,:)),'k','LineWidth',3)
hold on
for j = 1:(nC-1)
plot(omegaLoop,real(firstDSMode(j+1,:)),'LineWidth',2,'Color',c(j,:))
end
axis([0,omegaLoop(end),0,.014])
hold off
xlabel('$\omega$')
ylabel('$W_-(0)$')
legend(legendInfo,'Interpreter','latex','Location','northeast')

cleanfigure
matlab2tikz([imageFolder,'spod1-glegg','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

%% Plot second sound power output downstream mode
figure(3)
plot(omegaLoop,real(secondDSMode(1,:)),'k','LineWidth',3)
    legendInfo{1} = ['$C_{II} = ' num2str(0) '$']; 
hold on
for j = 1:(nC-1)
plot(omegaLoop,real(secondDSMode(j+1,:)),'LineWidth',2,'Color',c(j,:))
    legendInfo{j+1} = ['$C_{II} = ' num2str(C(j+1)) '$']; 
end
axis([0,omegaLoop(end),0,.014])
hold off
xlabel('$\omega$')
ylabel('$W_-(1)$')
legend(legendInfo,'Interpreter','latex','Location','northeast')

cleanfigure
matlab2tikz([imageFolder,'spod2-glegg','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

%% Plot difference in sound power output first mode
firstDSModeDel = 10*log10(SPOD(:,:,ceil(end/2))./SPOD(1,:,ceil(end/2)));
cutLocs1 = find(abs(SPOD(1,:,ceil(end/2)))>1e-5,1,'first');
firstDSModeDel(:,1:cutLocs1-1,ceil(end/2))=0;
secondDSModeDel = 10*log10(SPOD(:,:,ceil(end/2)+1)./SPOD(1,:,ceil(end/2)+1));
cutLocs2 = find(abs(SPOD(1,:,ceil(end/2)+1))>1e-5,1,'first');
secondDSModeDel(:,1:cutLocs2-1,ceil(end/2))=0;

figure(4)
clear legendInfo
for j = 1:(nC-1)
    plot(omegaLoop,real(firstDSModeDel(j+1,:)),'LineWidth',2,'Color',c(j,:))
    hold on
    legendInfo{j} = ['$C_{II} = ' num2str(C(j+1)) '$']; 
end

hold off
xlim([0,30])

xlabel('$\omega$');
ylabel('$L_{W_-}(0)$');
legend(legendInfo,'Interpreter','latex','Location','southwest')
cleanfigure
matlab2tikz([imageFolder,'spld1-glegg','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

%% Plot difference in sound power output second mode

figure(5)
clear legendInfo
for j = 1:(nC-1)
    plot(omegaLoop,real(secondDSModeDel(j+1,:)),'LineWidth',2,'Color',c(j,:))
    hold on
    legendInfo{j} = ['$C_{II} = ' num2str(C(j+1)) '$']; 
end

hold off
xlim([0,30])

xlabel('$\omega$');
ylabel('$L_{W_-}(1)$');
legend(legendInfo,'Interpreter','latex','Location','southwest')

cleanfigure
matlab2tikz([imageFolder,'spld2-glegg','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
