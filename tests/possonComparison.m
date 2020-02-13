vaneDist = .6;
stagAng = 40;

C1 = [.1,.1,3];
addpath(genpath('../'));
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


%l
%%
ADData=struct('spac',   [s,d],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       M,...                         %Mach number 
              'case',    'II', ...                       %Case
              'C',       0);                            %Coefficient of case                       
             
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
out = computeModes(newADData,newAAData,Modes);
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
figure(1)

plotField(h,newADData,newAAData,plotData,type)
%caxis([-30,30])