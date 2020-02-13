
addpath(genpath('../'));
imageFolder = '../../images/';

vaneDistDim = .6;
stagAngDim = 40;

dDim = vaneDistDim*sin(stagAngDim*pi/180);
sDim = vaneDistDim*cos(stagAngDim*pi/180);

omegaDim = 25*(1 + 1e-5i);
%l
%%
ADData=struct('spacDim',     [sDim,dDim],... %Dimensional blade spacing
              'chordDim',    1,...           %Blade length
              'c0',          340,...         %Speed of sound
              'M',           M,...           %Mach number 
              'Wdim',        0,...           %Spanwise background velocity
              'case',    3, ...              %Case
              'C',       1);                 %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omegaDim',    omegaDim,...                 %Time Frequency
               'kyDim',       [],...                       %Normal frequency
               'kxDim',       'gust',...                   %Tangential frequency
               'kzDim',       10,...                       %Spanwise frequency
               'sigmao',      3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Amp',         [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
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
x = [linspace(-1,0,100),xSurf,linspace(1,2,100)];
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

plotFieldScattered(h,newADData,newAAData,plotData,type)
%caxis([-30,30])