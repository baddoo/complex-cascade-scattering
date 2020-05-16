
addpath(genpath('../'));
imageFolder = '../../images/';

chordDim = 1;
stagAngDim = 30;
vaneDistDim = chordDim/cos(stagAngDim*pi/180);

dDim = vaneDistDim*sin(stagAngDim*pi/180);
sDim = vaneDistDim*cos(stagAngDim*pi/180);

%%
ADData=struct('spacDim',     [sDim,dDim],... %Dimensional blade spacing
              'chordDim',    chordDim,...           %Blade length
              'c0',          340,...         %Speed of sound
              'M',           0.2,...           %Mach number 
              'Wdim',        0,...                %Spanwise background velocity
              'case',        2, ...               %Case
              'C',           1*(-.2623+0.0244i)); %Coefficient of case                       
             % Issue is when C has a negative real part... Roots are found
             % fine, but associating which ones are in the upper or lower
             % half plane is difficult.
%% Aeroacoustic Data
AAData=struct( 'omegaDim',    [],...                 %Time Frequency
               'omega',       20 + 1e-3i,...                 %Time Frequency
               'kxDim',       [],...                   %Tangential frequency
               'kx',          20,...                   %Tangential frequency
               'ky',          [],...
               'kyDim',       [],...                       %Normal frequency
               'kzDim',       [],...                       %Spanwise frequency
               'kz',          0,...
               'Sigmao',      1*pi/4,...                  %Interblade phase angle in (x,y)-space
               'Sigma',       [],...
               'Amp',         [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',500,...                     %Truncation of kernel modes
             'dmodes',100,...                      %Number of duct mode
             'amodes',100);

[newADData,newAAData] = prepareData(ADData,AAData);      
tic
out = computeModes(newADData,newAAData,Modes);
toc
%%
data=computeCoefficients(newADData,newAAData,out);
Z = linspace(-4,6,200) + linspace(0,newADData.spac(3)).'*1i*exp(-1i*newADData.chie);
X = real(Z); Y = imag(Z);
plotData = struct('X',X,...
                  'Y',Y,...
                  'axisLimits',chordDim*[-3,3,-2,2]);
              
type = 'pressure';
phi = computeField(data,type);

figure(3)
plotFieldScattered(phi,newADData,newAAData,plotData,type)
caxis(.05*[-5,5])

%% Check BC
xgrid = 1 + sin(pi/2*(linspace(-1,1,100))); xgrid(1) = []; xgrid(end) = [];
vel = computeField(data,'vvelocity');

LW = 'LineWidth'; FS = 'FontSize';
upperPressure = phi(xgrid);
lowerPressure = phi(xgrid + newADData.spac(2) + 1i*newADData.spac(1)).*exp(-1i*newAAData.Sigmao);
pressureJump = upperPressure - lowerPressure;

incVVel = newAAData.Amp(2)*exp(1i*newAAData.kx/newADData.Beta^2*xgrid).*exp(-1i*newAAData.omega*newADData.M.^2/newADData.Beta.^2.*xgrid);
totalVVel = vel(xgrid) + incVVel;
figure(4)
plot(xgrid,real(newADData.C*pressureJump),LW,3)
hold on
plot(xgrid,real(2*totalVVel),LW,3)
hold off

legend('$C_{II}\times$ pressure', '$2 \times$ vertical velocity','Interpreter','Latex',FS,15)
