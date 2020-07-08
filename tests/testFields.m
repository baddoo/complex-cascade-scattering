%addpath(genpath('../'));
%imageFolder = '../../images/';

chordDim = 1;
stagAngDim = 20;
vaneDistDim = 0.8;

dDim = vaneDistDim*sin(stagAngDim*pi/180);
sDim = vaneDistDim*cos(stagAngDim*pi/180);

%%
ADData=struct('spacDim',     [sDim,dDim],... %Dimensional blade spacing
              'chordDim',    chordDim,...           %Blade length
              'c0',          340,...         %Speed of sound
              'M',           0.4,...           %Mach number 
              'Wdim',        0,...                %Spanwise background velocity
              'case',        2, ...               %Case
              'C',           -1); %Coefficient of case                       
%% Aeroacoustic Data
AAData=struct( 'omegaDim',    [],...                 %Time Frequency
               'omega',       16*(1 + 1e-5i),...                 %Time Frequency
               'kxDim',       [],...                   %Tangential frequency
               'kx',          6,...                   %Tangential frequency
               'ky',          [],...
               'kyDim',       [],...                       %Normal frequency
               'kzDim',       [],...                       %Spanwise frequency
               'kz',          0,...                         
               'Sigmao',      1*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Sigma',       [],...
               'Amp',         [nan,1,nan]);          %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',500,...                     %Truncation of kernel modes
             'dmodes',70,...                      %Number of duct mode
             'amodes',70);                      % Number of acoustic modes

[newADData,newAAData] = prepareData(ADData,AAData);      
tic
out = computeModes(newADData,newAAData,Modes);
disp(['It took ' num2str(toc) ' seconds to compute the modes.'])
%%
data=computeCoefficients(newADData,newAAData,out);
Z = linspace(-4,6,200) + linspace(0,newADData.spac(3)).'*1i*exp(-1i*newADData.chie);
X = real(Z); Y = imag(Z);
plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',chordDim*[-3,3,-2,2]);
              
type = 'pressure';
phi = @(z) computeField(z,data,type);

figure(3)
tic
plotFieldScattered(phi,newADData,newAAData,plotData,type)
disp(['It took ' num2str(toc) ' seconds to plot the solution over ' num2str(numel(Z)) ' points per period window.'])
caxis(.05*[-5,5])

%% Check BC
xgrid = 1 + sin(pi/2*(linspace(-1,1,100))); xgrid(1) = []; xgrid(end) = [];
vel = @(z) computeField(z,data,'vvelocity');
LW = 'LineWidth'; FS = 'FontSize';
upperPressure = phi(xgrid);
lowerPressure = phi(xgrid + newADData.spac(2) + 1i*newADData.spac(1)).*exp(-1i*newAAData.Sigmao);
pressureJump = upperPressure - lowerPressure;

incVVel = newAAData.Amp(2)*exp(1i*newAAData.kx/newADData.Beta^2*xgrid).*exp(-1i*newAAData.omega*newADData.M.^2/newADData.Beta.^2.*xgrid);
pot = @(z) computeField(z,data,'potential');
dpot = pot(xgrid) - pot(xgrid + newADData.spac(2) + 1i*newADData.spac(1)).*exp(-1i*newAAData.Sigma);
totalVVel = vel(xgrid) + incVVel;
figure(4)
plot(xgrid,abs(newADData.C*pressureJump),LW,3)
hold on
plot(xgrid,abs(2*totalVVel),LW,3)

hold off

legend('$C_{II}\times$ pressure', '$2 \times$ vertical velocity','Interpreter','Latex',FS,15)
