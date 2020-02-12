vaneDist = 0.6;
stagAng = 40;

imageFolder = '../../images/';

kVals = [40,5];
nk = numel(kVals);
mu2 = [0,-.01,-.02,-.03];
nMu = numel(mu2);

for m = 1:nMu
for l = 1:nk
%l
%%
ADData=struct('spac',[vaneDist*cos(stagAng*pi/180),vaneDist*sin(stagAng*pi/180)],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       .3,...                         %Mach number 
              'mu',       1*[.001,0,mu2(m)]);                           
              
%% Aeroacoustic Data
AAData=struct( 'omega',    [],...                      %Tangential Frequency
               'kn',       [],...                       %Normal frequency
               'k',        2*kVals(l)*(1+1e-3i),...
               'kz',       0,...                       %Spanwise frequency
               'sigmao',   3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1e-2*ADData.M.*ADData.c0,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',200,...                     %Truncation of kernel modes
             'dmodes',15,...                      %Number of duct mode
             'amodes',15);

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

figure(5)
plotField(h,newADData,newAAData,plotData,type)
he = ADData.spac(1);
axis([-4*he,12*he,0,7*he])
if l==1 
caxis(5e-1*[-1,1])
else 
caxis(1*[-1,1])    
end

latexPNG(['glegg-pres-',num2str(l),'-',num2str(m)],imageFolder,gca);
end
end
return
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