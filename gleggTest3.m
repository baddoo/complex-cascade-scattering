vaneDist = .6;
stagAng = 40;

nFreq = 10;
freq = linspace(.1,100,nFreq);
nFreq = 1;

addpath('matlab2tikz/src','cfe')
imageFolder = '../../images/';

%l
%%
ADData=struct('spac',[vaneDist*cos(stagAng*pi/180),vaneDist*sin(stagAng*pi/180)],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       .3,...                         %Mach number 
              'case',    'I', ...                       %Case
              'C',       .01);                            %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omega',    23.5*2*(1+1e-3i),...                      %Tangential Frequency
               'kn',       [],...                       %Normal frequency
               'k',        13.5*2*(1+1e-3i),...
               'kz',       0,...                       %Spanwise frequency
               'sigmao',   3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',200,...                     %Truncation of kernel modes
             'dmodes',15,...                      %Number of duct mode
             'amodes',15);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);  
data=computeCoefficients(newADData,newAAData,out);

%%
xSurf = (1+sin(pi/2*linspace(-1,1)))/2;

x = [linspace(-6,0),xSurf,linspace(1,7)];
y = newADData.spac(1)*(1+sin(linspace(-1,1,100)*pi/2))/2;
[X,Y] = meshgrid(x,y);

plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',[-3,3,-2,2]);              
            
newData=computeExponents(data,plotData);

type = 'pressure';

h=computeField(newData,type);
%%
figure(5)
%clf
plotField(h,newADData,newAAData,plotData,type)
caxis([-5,5])
