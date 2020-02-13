vaneDist = .8;
stagAng = 50;

addpath('matlab2tikz/src','cfe')

d = vaneDist*sin(stagAng*pi/180);
s = vaneDist*cos(stagAng*pi/180);

omega = [4,10,30,60]*(1+1e-3i);
%omega = 4*(1+1e-3i);
ky = -.1;
zy = -omega*ky;
M = 0.2;
m = 0;
kz = 0;
Beta = sqrt(1-M^2); delta = 1/Beta^2;
w=mysqrt(M*delta,kz/Beta); 

for l1 = [1,2,3,4]
%l
%%
ADData=struct('spac',   [s,d],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       M,...                         %Mach number 
              'case',    'II', ...                       %Case
              'C',       1e-2);                            %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omega',    omega(l1),...                      %Tangential Frequency
               'kn',       [],...                       %Normal frequency
               'k',        omega(l1),...
               'kz',       kz,...                       %Spanwise frequency
               'sigmao',   3*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',200,...                     %Truncation of kernel modes
             'dmodes',35,...                      %Number of duct mode
             'amodes',35);

[newADData,newAAData] = prepareData(ADData,AAData);         
out = computeModes(newADData,newAAData,Modes);        
%%
xSurf = (1+sin(pi/2*linspace(-1,1)))/2;

data=computeCoefficients(newADData,newAAData,out);
x = [linspace(-6,0,200),xSurf,linspace(1,7,200)];
%x = linspace(-6,8+0*newADData.spac(2),200);
y = newADData.spac(1)*(1+sin(linspace(-1,1,200)*pi/2))/2;
[X,Y] = meshgrid(x,y);

plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',[-3,3,-2.5,2.5]);              
            
newData=computeExponents(data,plotData);

type = 'pressure';

h=computeField(newData,type);
%%
omeg = linspace(0,2*pi,4+1);omeg(1) = [];

for l = 1:numel(omeg)
figure(6)

plotFieldTimeScattered(h,omeg(l),newADData,newAAData,plotData,type)
if l1 == 1; caxis([-.01,.01]); elseif l2 ==2; caxis([-.1,.1]); end
set(gca,'Position',[0 0 1 1]);
print(['../../../general-presentations/images/porous/porCase',num2str(l1),'plot',num2str(l)],'-dpng','-r50');
end

end