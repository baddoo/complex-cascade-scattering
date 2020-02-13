vaneDist = .6;
stagAng = 30;

C1 = [0.01,.5,3];
addpath('matlab2tikz/src','cfe')
%imageFolder = '../../images/';

d = vaneDist*sin(stagAng*pi/180);
s = vaneDist*cos(stagAng*pi/180);

omega = 40*(1+1e-3i);
ky = .1;
zy = -omega*ky;
M = 0.3;
m = 0;
kz = 0;
Beta = sqrt(1-M^2); delta = 1/Beta^2;
w=mysqrt(M*delta,kz/Beta); 

kx = 1/delta*mysqrt(omega*w,(s*Beta*zy+2*m*pi)/(s*Beta));

%if imag(kx)>1e-2, warning('The mode is not cut-on.'); return; end

sigma = kx*delta*d + ky*omega*(s*Beta);
sigmaTest = -(s*Beta)*zy + d*mysqrt(omega*w,(Beta*s*zy+2*m*pi)/(Beta*s));
sigmao = sigma - d*omega*M^2*delta;

for l1 = 1;%[1,2,3]
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
             'trunc',200,...                     %Truncation of kernel modes
             'dmodes',35,...                      %Number of duct mode
             'amodes',35);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);        
%%
xSurf = (1+sin(pi/2*linspace(-1,1)))/2;

data=computeCoefficients(newADData,newAAData,out);
x = [linspace(-6,0,200),xSurf,linspace(1,7,200)];
%x = linspace(-6,8+0*newADData.spac(2),200);
y = newADData.spac(1)*(1+sin(linspace(-1,1,200)*pi/2))/2;
[X,Y] = meshgrid(x,y);

plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',[-3,3,-2.5,2.5]);              
            
newData=computeExponents(data,plotData);

type = 'pressure';

h=computeField(newData,type);
%%
omeg = linspace(0,2*pi,40+1);omeg(1) = [];
for l = 1:numel(omeg)
figure(5)

plotFieldScattered(h.*exp(-1i*omeg(l)),newADData,newAAData,plotData,type)
caxis([-40,40])
set(gca,'Position',[0 0 1 1]);
print(['../../../general-presentations/images/animated-images/scatAnim',num2str(l)],'-dpng','-r50');
end

for l = 1:numel(omeg)
figure(6)

plotFieldTime(h,omeg(l),newADData,newAAData,plotData,type)
caxis([-40,40])
set(gca,'Position',[0 0 1 1]);
print(['../../../general-presentations/images/animated-images/totalAnim',num2str(l)],'-dpng','-r50');
end
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
%s = sqrt(newADData.spac(1).^2+newADData.spac(2)^2);
se = newADData.spac(3);

R = pi*abs(ZPa.*DLP./SQRTa./se);
T = pi*abs(ZMa.*DLM./SQRTa./se.*exp(-1i*LMa));
T(ceil(end/2))=T(ceil(end/2))+1;

modeNum = 2;

locs = find(abs(imag(LMa(floor(end/2)-modeNum+1:floor(end/2)+modeNum+1)))>1e-1);

num = numel(LMa);
x = -modeNum:modeNum;

figure(2)
%clf
b=bar(x,abs([R(floor(end/2)-modeNum+1:floor(end/2)+modeNum+1), ...
         T(floor(end/2)-modeNum+1:floor(end/2)+modeNum+1)]) ,...
         'FaceColor','flat');

cols = lines;   
     
b(1).CData(locs,:) = b(1).CData(locs,:)+.4;   
b(2).CData(locs,:) = b(2).CData(locs,:)+.4;     
xlabel('Mode number, m ')
ylabel('$\left| R_m \right|, \left| T_m \right|$')

cleanfigure
ylim([0,1.5])
%%matlab2tikz([imageFolder,'tr-',num2str(l1),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right');
%latexPNGNAR(['tr-',num2str(l1)],imageFolder,gca);
%close
end