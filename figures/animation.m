vaneDist = .6;
stagAng = 40;

addpath('matlab2tikz/src','cfe')
imageFolder = '../../images/';

d = vaneDist*sin(stagAng*pi/180);
s = vaneDist*cos(stagAng*pi/180);

omega = 100*(1 + 1e-5i);
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
C1 = .01;
%l
%%
ADData=struct('spac',   [s,d],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       M,...                         %Mach number 
              'case',    'II', ...                       %Case
              'C',       C1);                            %Coefficient of case                       
             
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
             'dmodes',15,...                      %Number of duct mode
             'amodes',15);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);

%%
xSurf = (1+sin(pi/2*linspace(-1,1)))/2;

data=computeCoefficients(newADData,newAAData,out);
x = [linspace(-6,0,300),xSurf,linspace(1,7,300)];
y = newADData.spac(1)*(1+sin(linspace(-1,1,300)*pi/2))/2;
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
omeg = linspace(0,2*pi,60+1);omeg(1) = [];
types = ["pressure","hvelocity"];
for j = 1:2
    h=computeField(newData,types(j));
for l = 1:numel(omeg)
figure(j)
plotFieldTime(h,omeg(l),newADData,newAAData,plotData,types(j))
if j ==1; caxis([-100,100]); elseif j ==2 ; caxis([-30,30]); end
set(gca,'Position',[0 0 1 1]);
print(strcat('animations/totalAnim',types(j),num2str(l,'%05.3d')),'-dpng','-r50');
l
end
end
return