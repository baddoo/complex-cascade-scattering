
addpath(genpath('../'));
imageFolder = '../../images/';

chordDim = 1;
vaneDistDim = chordDim/cos(stagAngDim*pi/180);
stagAngDim = 40;

dDim = vaneDistDim*sin(stagAngDim*pi/180);
sDim = vaneDistDim*cos(stagAngDim*pi/180);

%l
%%
ADData=struct('spacDim',     [sDim,dDim],... %Dimensional blade spacing
              'chordDim',    chordDim,...           %Blade length
              'c0',          340,...         %Speed of sound
              'M',           0.2,...           %Mach number 
              'Wdim',        0,...           %Spanwise background velocity
              'case',        2, ...              %Case
              'C',           -.2623-0.0244i);                 %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omegaDim',    [],...                 %Time Frequency
               'omega',       50*(1+1e-3i),...                 %Time Frequency
               'kxDim',       [],...                   %Tangential frequency
               'kx',          20,...                   %Tangential frequency
               'ky',          [],...
               'kyDim',       [],...                       %Normal frequency
               'kzDim',       [],...                       %Spanwise frequency
               'kz',          0,...
               'sigmao',      3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'sigma',       [],...
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
Z = linspace(-2,3,200) + linspace(0,newADData.spac(3)).'*1i*exp(-1i*newADData.chie);
X = real(Z); Y = imag(Z);
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
caxis(1*[-5,5])