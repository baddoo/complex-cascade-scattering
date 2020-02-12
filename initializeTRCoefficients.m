vaneDist = 0;
stagAng = 0;

nFreq = 10;
freq = linspace(.1,100,nFreq);
nFreq = 1;

%l
%%
ADData=struct('spac',[1/1.25+0*vaneDist*cos(stagAng*pi/180),0*tan(20*pi/180)+0*vaneDist*sin(stagAng*pi/180)],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       .2,...                         %Mach number 
              'C',       1,...
              'case',    'I');                           
              
%% Aeroacoustic Data
AAData=struct( 'omega',    2,...                      %Tangential Frequency
               'kn',       [],...                       %Normal frequency
               'k',        10*(1+0e-0i),...
               'kz',       0,...                       %Spanwise frequency
               'sigmao',   3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1e-2*ADData.M.*ADData.c0,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',100,...                     %Truncation of kernel modes
             'dmodes',15,...                      %Number of duct mode
             'amodes',15);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);        
%%
data=computeCoefficients(newADData,newAAData,out);

LPa = out.LPa; LMa = out.LMa; ZPa = out.ZPa; ZMa = out.ZMa; SQRTa = out.SQRTa;
data.comb = [1 0 1 0];
DLP = permute(D(permute(LPa,[3,2,1]),data),[3,2,1]);
data.comb = [0 1 0 1];
DLM = permute(D(permute(LMa,[3,2,1]),data),[3,2,1]);

Beta = newADData.Beta;
s = sqrt(newADData.spac(1).^2+newADData.spac(2)^2);
se = newADData.spac(3);

R3 = pi*abs(ZPa.*DLP./SQRTa)./se;
T3 = pi*abs(ZMa.*DLM./SQRTa./se.*exp(-1i*LMa));

R = permute(R3,[3,2,1]);
T = permute(T3,[3,2,1]);

figure(1)
bar(abs(R))

figure(2)
bar(abs(T))
 data=computeCoefficients(newADData,newAAData,out);
 x = [linspace(-6,0),xSurf,linspace(1,7)];
% %x = linspace(-6,8+0*newADData.spac(2),200);
 y = newADData.spac(1)*(1+sin(linspace(-1,1,100)*pi/2))/2;
 [X,Y] = meshgrid(x,y);
% 
 plotData = struct('X',X,...
                   'Y',Y,...
                  'axisLimits',[-3,3,-2,2]);              
%             
 newData=computeExponents(data,plotData);
% 
 type = 'pressure';
% 
 h=computeField(newData,type);

% figure(5)
plotField(h,newADData,newAAData,plotData,type)

%latexPNG(['posson-fig-11-',num2str(l)],imageFolder,gca);
%end
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