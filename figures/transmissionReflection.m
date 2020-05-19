
imageFolder = '../../../LaTeX/porous-jfm-r2s/images/';

chordDim = 1;
semiChordDim = chordDim/2;
vaneDistDim = 0.8*chordDim;
stagAngDim = 30;

dDim = vaneDistDim*sin(stagAngDim*pi/180);
sDim = vaneDistDim*cos(stagAngDim*pi/180);

c0 = 340;
M = 0.1;
Beta = sqrt(1-M^2);

% Define Rayleigh conductivity
RDim = .02*semiChordDim; % Define radius of aperture
howeReducedFreq = 2; % Set Howe's reduced frequency
omega = howeReducedFreq*semiChordDim/RDim; % Define my omega
%KR = 2*RDim*(0.096 - 1i*1.030); % 1
KR = 2*RDim*(1.252 - 1i*0.705);% 2 
%KR = 2*RDim*(-.243 - 1i*.264); % .5
%KR = 2*RDim*(1.015+1i*0.127); % 3.5
%KR = 0.928 + 1i*.0385; % Define Rayleigh conductivity from Howe 1996 tables
alphaH = [.01,.05]; % Define open area fractions

N = alphaH/(pi*RDim^2);

tildeKR = N*KR;
CII = 2*tildeKR/(1i*omega);
CII = [1e-4,CII];

d = dDim/semiChordDim;
s = sDim/semiChordDim*Beta;
del = sqrt(s^2+d^2); chie = atan(d/s);

Beta = sqrt(1-M^2);
kz = 0;
w=mysqrt(M/Beta^2,kz/Beta); % Need to spanwise flow terms.
Sigma = -pi/4;
fm = Sigma/del;
kx = -(-fm*sin(chie) - cos(chie).*mysqrt(omega*w,Sigma/del))*Beta^2;

for l1 = 1:3
%l
%%
ADData=struct('spacDim', [sDim,dDim],...%Blade Spacing [h,d]
              'chordDim',chordDim,...                          %Blade length
              'c0',      c0,...                        %Speed of sound
              'M',       M,...                         %Mach number 
              'Wdim',    0,...
              'case',    2, ...                       %Case
              'C',       CII(l1));                            %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omegaDim', [],...
               'omega',    omega*(1 + 1e-5i),...                      %Tangential Frequency
               'kx',       kx,...
               'kxDim',    [],...
               'ky',       [],...                       %Normal frequency
               'kyDim',    [],...
               'kz',       0,...                       %Spanwise frequency
               'kzDim',    [],...                       %Spanwise frequency
               'Sigma',    Sigma,...                  %Interblade phase angle in (x,y)-space
               'Sigmao',   [],...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',500,...                     %Truncation of kernel modes
             'dmodes',15,...                      %Number of duct mode
             'amodes',15);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);

%%
xSurf = (1+sin(pi/2*linspace(-1,1)))/2;

data=computeCoefficients(newADData,newAAData,out);

Z = linspace(-4,6,300) + linspace(0,newADData.spac(3),300).'*1i*exp(-1i*newADData.chie);
X = real(Z); Y = imag(Z);

plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',[-3,3,-2,2]);              
tic            
pres = computeField(data,'pressure');
toc


%%
figure(5)

plotField(pres,newADData,newAAData,plotData,'pressure')
caxis(1.2*[-100,100])
latexPNG(['trans-',num2str(l1)],imageFolder,gca);
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
T(ceil(end/2))=1+T(ceil(end/2));

modeNum = 2;

locs = find(abs(imag(LMa(floor(end/2)-modeNum+1:floor(end/2)+modeNum+1)))>1e-4);

num = numel(LMa);
x = -modeNum:modeNum;

 figure(2)
% Find modes that are cut off
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
ylim([0,2])
latexPNGNAR(['tr-',num2str(l1)],imageFolder,gca);
%close
end
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
ylim([0,1])
