
addpath(genpath('../'));
imageFolder = '../../images/';

chordDim = 1;
stagAngDim = 20;
vaneDistDim = chordDim/cos(stagAngDim*pi/180);

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
              'C',           (-.2623-0.0244i));                 %Coefficient of case                       
             % Issue is when C has a negative real part... Roots are found
             % fine, but associating which ones are in the upper or lower
             % half plane is difficult.
%% Aeroacoustic Data
AAData=struct( 'omegaDim',    [],...                 %Time Frequency
               'omega',       20 + 1e-3i,...                 %Time Frequency
               'kxDim',       [],...                   %Tangential frequency
               'kx',          15,...                   %Tangential frequency
               'ky',          [],...
               'kyDim',       [],...                       %Normal frequency
               'kzDim',       [],...                       %Spanwise frequency
               'kz',          0,...
               'Sigmao',      0*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Sigma',       [],...
               'Amp',         [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',500,...                     %Truncation of kernel modes
             'dmodes',35,...                      %Number of duct mode
             'amodes',35);

[newADData,newAAData] = prepareData(ADData,AAData);      
tic
out = computeModes(newADData,newAAData,Modes);
toc
%%
xSurf = chordDim*(1+sin(pi/2*linspace(-1,1)))/2; %xSurf(1) = []; xSurf(end) = [];

data=computeCoefficients(newADData,newAAData,out);
Z = linspace(-4,6,200) + linspace(0,newADData.spac(3)).'*1i*exp(-1i*newADData.chie);
X = real(Z); Y = imag(Z);
plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',chordDim*[-3,3,-2,2]);              
tic            
%newData=computeExponents(data,plotData);
toc
type = 'vvelocity';

phi = computeField(data,type);

figure(3)

plotFieldScattered(phi,newADData,newAAData,plotData,type)
caxis(.5*[-5,5])