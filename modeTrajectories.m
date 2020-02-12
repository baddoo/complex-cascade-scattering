s = 1; d = .3;
M = 0.3;

nPts = 100;
nDM = 5;
nAng = 8;

TPfin = zeros(nDM,nPts,nAng);
TMfin = zeros(nDM,nPts,nAng);

ang = (0:(nAng-1))/nAng*pi*2;

for l0 = 1:nAng
    C1 = exp(1i*ang(l0)).*10.^(linspace(-2,2,nPts));
for l1 = 1:nPts
ADData=struct('spac',   [s,d],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       M,...                         %Mach number 
              'case',    'II', ...                       %Case
              'C',       C1(l1));                            %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omega',    5,...                      %Tangential Frequency
               'kn',       [],...                       %Normal frequency
               'k',        3,...
               'kz',       0,...                       %Spanwise frequency
               'sigmao',   2*pi/3,...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
%% Information about Modes            
Modes=struct('comb',[1,1,1,1],...
             'trunc',nDM+1,...                     %Truncation of kernel modes
             'dmodes',nDM,...                      %Number of duct mode
             'amodes',5);

[newADData,newAAData] = prepareData(ADData,AAData);    

[TP3,TM3,asympGuess] = findDuctModes(newADData,newAAData,Modes);

TPfin(:,l1,l0) = permute(TP3(1:nDM),[3,1,2]);
TMfin(:,l1,l0) = permute(TM3(1:nDM),[3,1,2]);

end
end
%%

rigidModesP = mysqrt(newAAData.w.*newAAData.omega5,(0:nDM).*pi/newADData.spac(1));
rigidModesM = -rigidModesP;

TP2 = TPfin(:,l1,1);

for l1 = 1:(nPts-1)
v1 = TPfin(:,l1,1);
v2 = TPfin(:,l1+1,1);

d = v1-v2.';
[~,locs] = min(abs(d),[],1);

v2new = v2(locs);

TP2 = [TP2,v2new];

end


for l0 = 1:nAng

plot(TPfin(:,:,l0).','.')
hold on
plot(TMfin(:,:,l0).','.')

end

plot(real(rigidModesP),imag(rigidModesP),'ko','MarkerFaceColor','k')
plot(real(rigidModesM),imag(rigidModesM),'ko','MarkerFaceColor','k')


hold off
