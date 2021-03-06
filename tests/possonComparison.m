
addpath(genpath('../'));
%imageFolder = '../../images/';

chordDim = 1;
vaneDistDim = chordDim/cos(stagAngDim*pi/180);
stagAngDim = 20;

dDim = vaneDistDim*sin(stagAngDim*pi/180);
sDim = vaneDistDim*cos(stagAngDim*pi/180);

%l
%%
ADData=struct('spacDim',     [sDim,dDim],... %Dimensional blade spacing
              'chordDim',    chordDim,...           %Blade length
              'c0',          340,...         %Speed of sound
              'M',           0.2,...           %Mach number 
              'Wdim',        0,...           %Spanwise background velocity
              'case',        1, ...              %Case
              'C',           1e-3);                 %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omegaDim',    [],...                 %Time Frequency
               'omega',       12.5,...                 %Time Frequency
               'kxDim',       [],...                   %Tangential frequency
               'kx',          'gust',...                   %Tangential frequency
               'kyDim',       [],...                       %Normal frequency
               'kzDim',       0,...                       %Spanwise frequency
               'sigmao',      3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Amp',         [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',500,...                     %Truncation of kernel modes
             'dmodes',35,...                      %Number of duct mode
             'amodes',35);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);
%%
xSurf = chordDim*(1+sin(pi/2*linspace(-1,1)))/2; %xSurf(1) = []; xSurf(end) = [];

data=computeCoefficients(newADData,newAAData,out);
x = chordDim*[linspace(-4,0,100),xSurf,linspace(chordDim,6,100)];
y = chordDim*newADData.spac(1)*(1+sin(linspace(-1,1,50)*pi/2))/2;
[X,Y] = meshgrid(x,y);

plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',chordDim*[-3,3,-2,2]);              
tic            
newData=computeExponents(data,plotData);
toc
type = 'pressure';

h=computeField(newData,type);
%%
figure(1)

plotFieldScattered(h,newADData,newAAData,plotData,type)
caxis(1*[-1,1])