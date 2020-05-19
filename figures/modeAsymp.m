imageFolder = '../../../LaTeX/porous-jfm-r2s/images/';

M = .3;
stagAngDim = 30;
vaneDistDim = .5;
nPts = 100;
nDM = 5;
nAng = 8;
om = [5,10];
TPfin = zeros(2*nDM,nPts,nAng);

xSol = zeros(2*nDM,1);
ang = (0:(nAng-1))/nAng*2*pi;
C1vec = exp(1i*ang).*[0,10.^(linspace(-4,3,nPts-1))].';
for m = 2
    
for l0 = 1:nAng
    for l1 = 1:nPts
        C1 = C1vec(l1,l0);
ADData=struct('spacDim',   [vaneDistDim*cos(stagAngDim*pi/180),vaneDistDim*sin(stagAngDim*pi/180)],...%Blade Spacing
              'chordDim',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       M,...                         %Mach number 
              'Wdim',    0,...
              'case',    2, ...                       %Case
              'C',       C1);                            %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omega',    om(m)*(1+1e-9i),...                      %Tangential Frequency
               'omegaDim',[],...
               'kn',       [],...                       %Normal frequency
               'kx',       4*(1+1e-9i),...
               'kz',       0,...                       %Spanwise frequency
               'Sigmao',   2*pi/3,...                  %Interblade phase angle in (x,y)-space
               'Sigma',[],'kxDim',[],'ky',[],'kyDim',[],'kzDim',[],...
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
% Prepare data
[newADData,newAAData] = prepareData(ADData,AAData);    

    mu = newADData.mu;
    dPG = newADData.spac(2); sPG = newADData.spac(1); omega = newAAData.omega;
    delPG = newADData.spac(3);
    w = newAAData.w;
    Sigma = newAAData.Sigma;
    fm = (Sigma - 2*pi*(-nDM:nDM))/delPG;

    if l1 == 1
    
    rDP = mysqrt(omega*w,(0:(nDM-1)).'.*pi/sPG);
    rDM = - rDP;
    rAP = -fm*dPG./delPG + sPG/delPG.*mysqrt(w*omega,fm);
    rAM = -fm*dPG./delPG - sPG/delPG.*mysqrt(w*omega,fm);
    rigidAModes = [rAP;rAM];
    rigidDModes = [rDP;rDM];
    TPfin(:,l1,l0) = rigidDModes;

    else
  
    f = @(xVar) 1./KNumlogD(xVar,newADData,newAAData);
    logD = @(xVar) KNumlogD(xVar,newADData,newAAData);
    tol = 1e-5;
    TPfin(:,l1,l0) = newtonKRoots(TPfin(:,l1-1,l0),f,logD,tol,[]);
    
    end
[l0,l1]
    end
end


%% Plots
figure(1)
clf
hold on

recCol1 = [255,240,240]/255;
recCol2 = [240,240,255]/255;
rectangle('Position',[-500,0,1000,500],'FaceColor',recCol1,'EdgeColor',recCol1)
rectangle('Position',[-500,-500,1000,500],'FaceColor',recCol2,'EdgeColor',recCol2)
probLoc = find(ang==pi*3/2);
loc1 = find(abs(TPfin(1,:,probLoc))>10);

if ~isempty(loc1); TPfin(1,loc1(1):end,probLoc)=nan; end

c = hsv(nAng);
for l0 = 1:nAng
plot(TPfin(:,:,l0).','Color',c(l0,:),'LineWidth',1)
end
plot(rigidDModes,'ko','MarkerFaceColor','k')%,'MarkerSize',8)
plot(rigidAModes,'ks','MarkerFaceColor','k')%,'MarkerSize',8)
% Convected mode
cMode = -omega/newADData.Beta.^2;
plot(real(cMode),imag(cMode),'k^','MarkerFaceColor','k')%,'MarkerSize',25)
%set(gca,'color',.9*[1,1,1]);
axis equal

hold off
axis([-15,10,-17,17])
xlabel('$\Re[\theta]$');%,'Interpreter','Latex')
ylabel('$\Im[\theta]$');%,'Interpreter','LaTeX')
set(gca,'Layer','top')

% cleanfigure
% matlab2tikz([imageFolder,'mode-traj-',num2str(m),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
%                  'extraaxisoptions','ylabel shift = -.4cm,,axis on top');
%%             
 if m==1
figure(2)

Csmall = -logspace(-4,0);
Clarge = -logspace(-1,3);
nAsym = (0:(nDM-1)).';
dN = ones(nDM,1); dN(1)=2;
zeta = @(zVar) mysqrt(omega*w,zVar);
smallRootsP = rDP - Csmall.*1i.*(rDP-cMode)./(dN.*sPG.*rDP).*(1-(-1).^nAsym.*cos(dPG.*rDP+Sigma));
smallRootsM = rDM - Csmall.*1i.*(rDM-cMode)./(dN.*sPG.*rDM).*(1-(-1).^nAsym.*cos(dPG.*rDM+Sigma));

zP = (dPG*rAP + delPG*fm)/sPG;
zM =-(dPG*rAM + delPG*fm)/sPG;

largeRootsP = rAP - 1./Clarge.'.*1i.*(zP.^2)./((rAP-cMode).*delPG.*zeta(fm));
largeRootsM = rAM + 1./Clarge.'.*1i.*(zM.^2)./((rAM-cMode).*delPG.*zeta(fm));
largeCmode = cMode - 1./Clarge*1i*zeta(cMode).*sin(sPG*zeta(cMode))./(cos(sPG*zeta(cMode))-cos(dPG*cMode+Sigma));
recCol1 = [255,240,240]/255;
recCol2 = [240,240,255]/255;

semilogx(abs(C1vec(:,end/2+1)),imag(TPfin(:,:,end/2+1)),'k','LineWidth',2)
hold on
semilogx(abs(Csmall),imag(smallRootsP).','b--')
semilogx(abs(Clarge),imag(largeCmode).','r--')
semilogx(abs(Clarge),imag(largeRootsP(:,[1,4,5,6,7,8])),'r--')
hold off
axis([1e-4,1e2,-.01,16])
xlabel('$C_{II}$');%,'Interpreter','Latex')
ylabel('$\Im[\theta]$');%,'Interpreter','LaTeX')
set(gca,'color',recCol1);

% cleanfigure
% matlab2tikz([imageFolder,'asympU.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
%                  'extraaxisoptions','axis on top');


figure(3)
semilogx(abs(C1vec(:,1)),imag(TPfin(:,:,end/2+1)),'k','LineWidth',2)
hold on
semilogx(abs(Csmall),imag(smallRootsM).','b--')
semilogx(abs(Clarge),imag(largeCmode).','r-.')
semilogx(abs(Clarge),imag(largeRootsM(:,[2,5,6,7,8,11])),'r--')
hold off
axis([1e-4,1e2,-16,0.01])
xlabel('$C_{II}$');%,'Interpreter','Latex')
ylabel('$\Im[\theta]$');%,'Interpreter','LaTeX')
set(gca,'color',recCol2);

% cleanfigure
% matlab2tikz([imageFolder,'asympL.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
%                  'extraaxisoptions','axis on top');

 end
end
