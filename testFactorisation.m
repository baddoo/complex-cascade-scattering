vaneDist = 0;
stagAng = 0;

nFreq = 10;
freq = linspace(.1,100,nFreq);
lift = zeros(size(freq));

nFreq = 1;

%for l = 1:nFreq;
%l
%%
ADData=struct('spac',[1+0*vaneDist*cos(stagAng*pi/180),tan(20*pi/180)+0*vaneDist*sin(stagAng*pi/180)],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       .2,...                         %Mach number 
              'mu',       1e-2);                           
              
%% Aeroacoustic Data
AAData=struct( 'omega',    [],...                      %Tangential Frequency
               'kn',       [],...                       %Normal frequency
               'k',        5/.2*(1+1e-9i),...
               'kz',       0,...                       %Spanwise frequency
               'sigmao',   3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1e-2*ADData.M.*ADData.c0,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',100,...                     %Truncation of kernel modes
             'dmodes',15,...                      %Number of duct mode
             'amodes',25);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);        
%%
xSurf = (1+sin(pi/2*linspace(-1,1)))/2;

data=computeCoefficients(newADData,newAAData,out);
x = [linspace(-6,0),xSurf,linspace(1,7)];
%x = linspace(-6,8+0*newADData.spac(2),200);
y = newADData.spac(1)*(1+sin(linspace(-1,1,100)*pi/2))/2;
[X,Y] = meshgrid(x,y);

plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',[-3,3,-2,2]);              
            
newData=computeExponents(data,plotData);


type = 'pressure';

h=computeField(newData,type);

figure(1)
plotField(h,newADData,newAAData,plotData,type)

%% Surface pressure

y = 0;
[X,Y] = meshgrid(xSurf,y);
surfPlotDataU = struct('X',X,...
                  'Y',Y);              
surfPlotDataL = struct('X',X+newADData.spac(2),...
                  'Y',Y+newADData.spac(1));        
              
surfUNewData=computeExponents(data,surfPlotDataU);
surfLNewData=computeExponents(data,surfPlotDataL);

hSurfU=computeField(surfUNewData,'pressure');
hSurfL=computeField(surfLNewData,'pressure');

figure(2)
plot(xSurf,real(hSurfU))
hold on
plot(xSurf,real(hSurfL.*exp(-1i*newAAData.sigmao)))
hold off
%lift(l) = abs(D(0,data));

%end
%factorizeKRadInt(newADData,newAAData,Modes);