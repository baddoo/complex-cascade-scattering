vaneDist = .6;
stagAng = 40;

C1 = [.1,.1,3];
addpath('matlab2tikz/src','cfe')
imageFolder = '../../images/';

d = vaneDist*sin(stagAng*pi/180);
s = vaneDist*cos(stagAng*pi/180);

omega = 25*(1 + 1e-15i);
ky = .1;
zy = -omega*ky;
M = 0.3;
m = 0;
kz = 0;
Beta = sqrt(1-M^2); delta = 1/Beta^2;
w=mysqrt(M*delta,kz/Beta); 

kx = 1/delta*mysqrt(omega*w,(s*Beta*zy+2*m*pi)/(s*Beta));

if imag(kx)>1e-2, warning('The mode is not cut-on.'); return; end

sigma = kx*delta*d + ky*omega*(s*Beta);
sigmaTest = -(s*Beta)*zy + d*mysqrt(omega*w,(Beta*s*zy+2*m*pi)/(Beta*s));
sigmao = sigma - d*omega*M^2*delta;

for l1 = 1%[1,2,3]
%l
%%
ADData=struct('spac',   [s,d],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       M,...                         %Mach number 
              'case',    'II', ...                       %Case
              'C',       C1(l1));                            %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omega',    omega,...                      %Tangential Frequency
               'kn',       [],...                       %Normal frequency
               'k',        kx,...
               'kz',       kz,...                       %Spanwise frequency
               'sigmao',   sigmao,...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',500,...                     %Truncation of kernel modes
             'dmodes',35,...                      %Number of duct mode
             'amodes',35);

[newADData,newAAData] = prepareData(ADData,AAData);         
%tic
%out = computeModes(newADData,newAAData,Modes); toc  
%tic
out = computeModes2(newADData,newAAData,Modes);
%toc
%%
xSurf = (1+sin(pi/2*linspace(-1,1)))/2;

data=computeCoefficients(newADData,newAAData,out);
x = [linspace(-6,0,100),xSurf,linspace(1,7,100)];
y = newADData.spac(1)*(1+sin(linspace(-1,1,50)*pi/2))/2;
[X,Y] = meshgrid(x,y);

plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',[-3,3,-2,2]);              
tic            
newData=computeExponents(data,plotData);
toc
type = 'pressure';

h=computeField(newData,type);
%%
figure(5)

plotField(h,newADData,newAAData,plotData,type)
%caxis([-30,30])
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
%s = sqrt(newADData.spac(1).^2+newADData.spac(2)^2);
se = newADData.spac(3);

R = pi*abs(ZPa.*DLP./SQRTa./se);
T = pi*abs(ZMa.*DLM./SQRTa./se.*exp(-1i*LMa));
T(ceil(end/2))=T(ceil(end/2))+1;

modeNum = 2;

locs = find(abs(imag(LMa(floor(end/2)-modeNum+1:floor(end/2)+modeNum+1)))>1e-1);

num = numel(LMa);
x = -modeNum:modeNum;

% figure(2)
% %clf
% b=bar(x,abs([R(floor(end/2)-modeNum+1:floor(end/2)+modeNum+1), ...
%          T(floor(end/2)-modeNum+1:floor(end/2)+modeNum+1)]) ,...
%          'FaceColor','flat');
% 
% cols = lines;   
%      
% b(1).CData(locs,:) = b(1).CData(locs,:)+.4;   
% b(2).CData(locs,:) = b(2).CData(locs,:)+.4;     
% xlabel('Mode number, m ')
% ylabel('$\left| R_m \right|, \left| T_m \right|$')
% 
% cleanfigure
% ylim([0,1.5])
% %%matlab2tikz([imageFolder,'tr-',num2str(l1),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
% latexPNGNAR(['tr-',num2str(l1)],imageFolder,gca);
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