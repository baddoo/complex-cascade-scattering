%addpath(genpath('../'));
imageFolder = [];

chordDim = 1;
semiChordDim = chordDim/2;
vaneDistDim = 0.6*chordDim;
stagAngDim = 30;

dDim = vaneDistDim*sin(stagAngDim*pi/180);
sDim = vaneDistDim*cos(stagAngDim*pi/180);

c0 = 340;
M = 0.2;
Beta = sqrt(1-M^2);

% Define Rayleigh conductivity
RDim = .02*semiChordDim; % Define radius of aperture
howeReducedFreq = .5; % Set Howe's reduced frequency
omega = howeReducedFreq*semiChordDim/RDim*(1+1e-3i); % Define my omega
%KR = 2*RDim*(0.096 - 1i*1.030); % 1
%KR = 2*RDim*(1.252 - 1i*0.705);% 2 
KR = 2*RDim*(-.243 - 1i*.264); % .5
%KR = 2*RDim*(1.015+1i*0.127); % 3.5
%KR = 0.928 + 1i*.0385; % Define Rayleigh conductivity from Howe 1996 tables
alphaH = [.05,.1]; % Define open area fractions

N = alphaH/(pi*RDim^2);

CII = -2*alphaH*KR/(1i*omega*pi*RDim^2);
CII = [1e-4,CII];

sigma = -3*pi/4;
kx = 5;

ADData=struct('spacDim', [sDim,dDim],...%Blade Spacing [h,d]
              'chordDim',chordDim,...                          %Blade length
              'c0',      c0,...                        %Speed of sound
              'M',       M,...                         %Mach number 
              'Wdim',    0,...
              'case',    2, ...                       %Case
              'C',       -.5i);                            %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omegaDim', [],...
               'omega',    5*(1 + 1e-7i),...                      %Tangential Frequency
               'kx',       5*(1 + 1e-7i),...
               'kxDim',    [],...
               'ky',       [],...                       %Normal frequency
               'kyDim',    [],...
               'kz',       0,...                       %Spanwise frequency
               'kzDim',    [],...                       %Spanwise frequency
               'sigma',    sigma,...                  %Interblade phase angle in (x,y)-space
               'sigmao',   [],...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',100,...                     %Truncation of kernel modes
             'dmodes',35,...                      %Number of duct mode
             'amodes',35);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);

%%
xSurf = (1+sin(pi/2*linspace(-1,1)))/2;

data=computeCoefficients(newADData,newAAData,out);

Z = linspace(-2,2+0*newADData.spac(2),100) + linspace(0,newADData.spac(3),100).'*1i*exp(-1i*newADData.chie);
%Z = linspace(0,newADData.spac(2))+1i*linspace(0,newADData.spac(1)-.01,500)';
X = real(Z); Y = imag(Z);

plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',[-3,3,-2,2]);              
tic            
newData=computeExponents(data,plotData);
toc
type = 'potential';

% Returns a function handle corresponding to the parameters and type.
h = @(zVar) evaluateField(zVar,newData,type);

pcolor(X,Y,real(h(Z))); shading interp
caxis(.1*[-1,1])
return
%%
figure(5)

% Returns a function 
plotFieldTest(h(Z),newADData,newAAData,plotData,type)
caxis(.005*[-100,100])

return

%%
LPa = permute(out.LPa,[3,2,1]);
LMa = permute(out.LMa,[3,2,1]);
ZPa = permute(out.ZPa,[3,2,1]);
ZMa = permute(out.ZMa,[3,2,1]); 
SQRTa = permute(out.SQRTa,[3,2,1]);

data.comb = [1 0 1 0];
DLP = D(LPa,data);
data.comb = [0 1 0 1];
DLM = D(LMa,data);

Beta = newADData.Beta;
se = newADData.spac(3);

R = pi*abs(ZPa.*DLP./SQRTa./se);
T = pi*abs(ZMa.*DLM./SQRTa./se.*exp(-2i*LMa));
T(ceil(end/2))=T(ceil(end/2));

modeNum = 2;

locs = find(abs(imag(LMa(floor(end/2)-modeNum+1:floor(end/2)+modeNum+1)))>1e-2);

num = numel(LMa);
x = -modeNum:modeNum;

 figure(2)
% %clf
 b=bar(x,abs([R(floor(end/2)-modeNum+1:floor(end/2)+modeNum+1), ...
          T(floor(end/2)-modeNum+1:floor(end/2)+modeNum+1)]) ,...
          'FaceColor','flat');
% 
 cols = lines;   
%      
 b(1).CData(locs,:) = b(1).CData(locs,:)+.4;   
 b(2).CData(locs,:) = b(2).CData(locs,:)+.4;     
 xlabel('Mode number, m ')
 ylabel('$\left| R_m \right|, \left| T_m \right|$','Interpreter','Latex')
% 
% cleanfigure
 ylim([0,2])
% %%matlab2tikz([imageFolder,'tr-',num2str(l1),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
latexPNGNAR(['tr-',num2str(l1)],imageFolder,gca);
%close

%%
nC = 200;
C1 = linspace(0.,100,nC/2+1); C1(1)=[];
C2 = 10.^linspace(-2,2,nC/2);
C = sort([C1,C2]);
%C(C==0)=[];
%C = 10.^linspace(-2,0,nC);
Ts = zeros(2,nC);
Rs = zeros(2,nC);
return
%%
for l = 1:nC
l
ADData.C = C(l);    
[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);        
data=computeCoefficients(newADData,newAAData,out);

LPa = permute(out.LPa,[3,2,1]);
LMa = permute(out.LMa,[3,2,1]);
ZPa = permute(out.ZPa,[3,2,1]);
ZMa = permute(out.ZMa,[3,2,1]); 
SQRTa = permute(out.SQRTa,[3,2,1]);

data.comb = [1 0 1 0];
DLP = D(LPa,data);
data.comb = [0 1 0 1];
DLM = D(LMa,data);

Beta = newADData.Beta;
%s = sqrt(newADData.spac(1).^2+newADData.spac(2)^2);
se = newADData.spac(3);

R = pi*abs(ZPa.*DLP./SQRTa./se);
T = pi*abs(ZMa.*DLM./SQRTa./se.*exp(-1i*LMa));
 
% Normalise
Rs(:,l) = R(ceil(end/2):ceil(end/2)+1);
Ts(:,l) = T(ceil(end/2):ceil(end/2)+1);
end

%%
blues = cols(1,:);
reds = cols(2,:);
figure(3)
loglog(abs(C),abs(Rs(1,:)./Rs(1,1)),'-','Color',blues,'LineWidth',3)
hold on
%loglog(abs(C),abs((Rs(1,:)-0*Rs(1,end-1))./Rs(1,1)),'--', 'Color',blues,'LineWidth',3)
loglog(abs(C),abs(Ts(1,:)./Ts(1,1)),'-','Color',reds,'LineWidth',3)
%loglog(abs(C),abs((Ts(1,:)-0*Ts(1,end-1))./Ts(1,1)),'--','Color',reds,'LineWidth',3)
xlabel('$C_{II}$')
ylabel('$\left| R_0 \right|, \left| T_0 - 1 \right|$')

hold off
ylim([0,1])
matlab2tikz([imageFolder,'trc.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

%%
figure(4)
semilogx(abs(C),abs(Rs(1,:)./Rs(1,1)),'-','Color',blues,'LineWidth',3)
hold on
%semilogx(abs(C),abs(Rs(1,:)./Rs(1,1)),'--', 'Color',blues,'LineWidth',3)
%semilogx(abs(C),abs(Ts(2,:))./Ts(2,1),'-','Color',reds,'LineWidth',3)
semilogx(abs(C),abs(Ts(1,:)./Ts(1,1)),'-','Color',reds,'LineWidth',3)
hold off

xlabel('$C_{II}$')
ylabel('$\left| R_0 \right|, \left| T_0-1 \right|$')
matlab2tikz([imageFolder,'trc-log.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');

%xlim([.0,.1])
ylim([0,1])
%%axis([-inf,1,1e-2,1])

function F = root(x,del,w,omega,chie,d,s)

f = (x(3)-2*pi)/del;
mysq = mysqrt(w*omega,f);

F(1) = (x(1) - (-f*sin(chie) - cos(chie)*mysq));
F(2) = (x(2) - (-(x(3) + d*x(1)-2*pi))/s);
F(3) = (x(3) - (-d*x(1) - s*x(2)));

end