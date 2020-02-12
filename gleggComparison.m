 vaneDist = 0.6;
stagAng = 40;

nFreq = 500;
C = [0.001,0.01,.1,1];
nC = numel(C);
kLoop = linspace(1,20,nFreq+1); kLoop(1)=[];
lift = zeros(nC,nFreq);
myAModes = 35;
SPOU = zeros(nC,nFreq,2*myAModes+1);
SPOD = zeros(nC,nFreq,2*myAModes+1);
addpath('matlab2tikz/src')
imageFolder = '../../images/';
for m = 1:nC
for l = 1:nFreq
disp(['Loop number ' num2str(l)])
%%
ADData=struct('spac',[vaneDist*cos(stagAng*pi/180),vaneDist*sin(stagAng*pi/180)],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       .3,...                         %Mach number 
              'case',    'II',...
              'C',       C(m));                           
              
%% Aeroacoustic Data
AAData=struct( 'omega',    2*kLoop(l)*(1+1e-3i),...                      %Tangential Frequency
               'kn',       [],...                       %Normal frequency
               'k',        2*kLoop(l)*(1+1e-3i),...
               'kz',       0,...                       %Spanwise frequency
               'sigmao',   3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',500,...                     %Truncation of kernel modes
             'dmodes',35,...                      %Number of duct mode
             'amodes',myAModes);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes2(newADData,newAAData,Modes);    
data=computeCoefficients(newADData,newAAData,out);

%%
lift(m,l) = 4i*kLoop(l).*D(-2*kLoop(l).*newADData.M^2.*newADData.delta,data)./newAAData.Amp5(2);
LPa = out.LPa; LMa = out.LMa; ZPa = out.ZPa; ZMa = out.ZMa; SQRTa = out.SQRTa;

data.comb = [1 0 1 0];
DLP = permute(D(permute(LPa,[3,2,1]),data),[3,2,1]);
data.comb = [0 1 0 1];
DLM = permute(D(permute(LMa,[3,2,1]),data),[3,2,1]);

Beta = newADData.Beta;
s = sqrt(ADData.spac(1).^2+ADData.spac(2)^2);
se = newADData.spac(3);

spoU = pi^2*4*Beta^2*kLoop(l).*abs(ZPa.*DLP).^2.*real(1./SQRTa)./s./abs(newAAData.Amp5(2).^2);
spoD = pi^2*4*Beta^2*kLoop(l).*abs(ZMa.*DLM).^2.*real(1./SQRTa)./s./abs(newAAData.Amp5(2).^2);
%spoD = abs(ZMa.^2);
%spoD = real((abs(ZMa.*DLM).^2)./SQRTa);

SPOU(m,l,:) = spoU;
SPOD(m,l,:) = spoD;

end
end
%% Plot lift
figure(1)
plot(kLoop,abs(lift(1,:)),'k','LineWidth',3);
hold on
plot(kLoop,abs(lift(2:end,:)),'LineWidth',2);
hold off
ylim([0,1.5])
xlim([0,20])


xlabel('$\omega$');%,'Interpreter','Latex')
ylabel('$\left|C_p\right|$');%,'Interpreter','LaTeX')
%newAAData.k = 2*kLoop;
%plotCutOnFrequencies(abs(lift(1,:,:).'),newADData,newAAData,Modes);
hold off
cleanfigure
matlab2tikz([imageFolder,'ul-glegg','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

%% Plot SPOD
firstDSMode = SPOD(:,:,ceil(end/2));
secondDSMode = SPOD(:,:,ceil(end/2)+1);
thirdDSMode = SPOD(:,:,ceil(end/2)+2);

figure(2)
plot(kLoop,real(firstDSMode(1,:)),'k','LineWidth',3)
hold on
plot(kLoop,real(firstDSMode(2:end,:)),'LineWidth',2)
h = gca;
set(h,'ColorOrderIndex',1);
plot(kLoop,real(secondDSMode(1,:)),'k--','LineWidth',3)
plot(kLoop,real(secondDSMode(2:end,:)),'--','LineWidth',2)
h = gca;
set(h,'ColorOrderIndex',1);
plot(kLoop,real(thirdDSMode(1,:)),'k-.','LineWidth',3)
plot(kLoop,real(thirdDSMode(2:end,:)),'-.','LineWidth',2)
axis([0,kLoop(end),0,.012])
hold off
xlabel('$\omega$')%,'Interpreter','Latex')
ylabel('Downstream sound power output')%,'Interpreter','LaTeX')
cleanfigure
matlab2tikz([imageFolder,'spod-glegg','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');


%% Plot SPOD
firstUSMode = SPOU(:,:,ceil(end/2));
secondUSMode = SPOU(:,:,ceil(end/2)+1);
thirdUSMode = SPOU(:,:,ceil(end/2)+2);

figure(3)
plot(kLoop,real(firstUSMode(1,:)),'k','LineWidth',3)
hold on
plot(kLoop,real(firstUSMode(2:end,:)),'LineWidth',2)
h = gca;
set(h,'ColorOrderIndex',1);
plot(kLoop,real(secondUSMode(1,:)),'k--','LineWidth',3)
plot(kLoop,real(secondUSMode(2:end,:)),'--','LineWidth',2)
h = gca;
set(h,'ColorOrderIndex',1);
plot(kLoop,real(thirdUSMode(1,:)),'k-.','LineWidth',3)
plot(kLoop,real(thirdUSMode(2:end,:)),'-.','LineWidth',2)
axis([0,kLoop(end),0,.025])
hold off

xlabel('$\omega$');%,'Interpreter','Latex')
ylabel('Upstream sound power output, $W$');%,'Interpreter','LaTeX')

cleanfigure
matlab2tikz([imageFolder,'spou-glegg','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

%% Plot SPOUDelta
firstUSModeDel = 10*log10(SPOU(:,:,ceil(end/2))./SPOU(1,:,ceil(end/2)));
firstUSModeDel(abs(SPOU(:,:,ceil(end/2)))<1e-10)=0;
secondUSModeDel = 10*log10(SPOU(:,:,ceil(end/2)+1)./SPOU(1,:,ceil(end/2)+1));
secondUSModeDel(abs(SPOU(:,:,ceil(end/2)+1))<1e-8)=0;
thirdUSModeDel = 10*log10(SPOU(:,:,ceil(end/2)+2)./SPOU(1,:,ceil(end/2)+2));
thirdUSModeDel(abs(SPOU(:,:,ceil(end/2)+2))<1e-8)=0;

figure(4)
plot(kLoop,real(firstUSModeDel(2:end,:)),'LineWidth',2)
hold on
h = gca;
set(h,'ColorOrderIndex',1);
plot(kLoop,real(secondUSModeDel(2:end,:)),'--','LineWidth',2)
h = gca;
set(h,'ColorOrderIndex',1);
plot(kLoop,real(thirdUSModeDel(2:end,:)),'-.','LineWidth',2)
%axis([0,kLoop(end),0,.025])
hold off
xlim([0,30])

xlabel('$\omega$');%,'Interpreter','Latex')
ylabel('Change in upstream sound power output');%,'Interpreter','LaTeX')

cleanfigure
matlab2tikz([imageFolder,'splu-glegg','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

%% Plot SPODDelta
%SPOD(abs(SPOD)<1e-5)=0;
%firstDSModeDel = 10*log10(sum(SPOD,3)./sum(SPOD(1,:,:),3));
firstDSModeDel = 10*log10(SPOD(:,:,ceil(end/2))./SPOD(1,:,ceil(end/2)));
firstDSModeDel(abs(SPOD(:,:,1+ceil(end/2)))<1e-12)=0;
secondDSModeDel = 10*log10(SPOD(:,:,ceil(end/2)+1)./SPOD(1,:,ceil(end/2)+1));
secondDSModeDel(abs(SPOD(:,:,ceil(end/2)+1))<1e-5)=0;
thirdDSModeDel = 10*log10(SPOD(:,:,ceil(end/2)+2)./SPOD(1,:,ceil(end/2)+2));
thirdDSModeDel(abs(SPOD(:,:,ceil(end/2)+2))<1e-5)=0;


figure(5)
%plot(kLoop,real(firstDSModeDel(1,:)),'k','LineWidth',3)
%hold on
plot(kLoop,real(firstDSModeDel(2:end,:)),'LineWidth',2)
hold on
h = gca;
set(h,'ColorOrderIndex',1);
plot(kLoop,real(secondDSModeDel(2:end,:)),'--','LineWidth',2)
h = gca;
set(h,'ColorOrderIndex',1);
plot(kLoop,real(thirdDSModeDel(2:end,:)),'-.','LineWidth',2)
%axis([0,kLoop(end),0,.025])
hold off
xlim([0,30])

xlabel('$\omega$');%,'Interpreter','Latex')
ylabel('Downstream sound power level, $L_W$');%,'Interpreter','LaTeX')

cleanfigure
matlab2tikz([imageFolder,'spld-glegg','.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
