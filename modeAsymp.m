imageFolder = '../../images/';


del = 1;
M = .3;
chi = [10,30,30]*pi/180;

nPts = 200;
nDM = 5;
nAng = 8;
om = [5,10];
TPfin = zeros(2*nDM,nPts,nAng);

xSol = zeros(2*nDM,1);
ang = (0:(nAng-1))/nAng*pi*2;
C1vec = exp(1i*ang).*[0,10.^(linspace(-4,3,nPts-1))].';
for m = 1:2
    
for l0 = 1:nAng
    for l1 = 1:nPts
        C1 = C1vec(l1,l0);
ADData=struct('spac',   [del*cos(chi(2)),del*sin(chi(2))],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       M,...                         %Mach number 
              'case',    'II', ...                       %Case
              'C',       C1);                            %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omega',    om(m)*(1+1e-9i),...                      %Tangential Frequency
               'kn',       [],...                       %Normal frequency
               'k',        4*(1+1e-9i),...
               'kz',       0,...                       %Spanwise frequency
               'sigmao',   2*pi/3,...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
% Prepare data
[newADData,newAAData] = prepareData(ADData,AAData);    

    mu = newADData.mu;
    dPG = newADData.spac(2); sPG = newADData.spac(1); omega = newAAData.omega5;
    delPG = newADData.spac(3);
    w = newAAData.w;
    sigma = newAAData.sigma5;
  F= @(x) mysqrt(omega*w,x).*sin(sPG*mysqrt(omega*w,x))./(cos(sPG*mysqrt(omega*w,x)) - cos(dPG*x + sigma)) ...
        + (-mu(3).*x.^2 -1i*mu(2).*x + mu(1));
fm = (sigma - 2*pi*(-nDM:nDM))/delPG;

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
%loc2 = find(abs(TPfin(2,:,probLoc))>10);
if ~isempty(loc1); TPfin(1,loc1(1):end,probLoc)=nan; end
%TPfin(2,loc2(2):end,probLoc)=nan;

c = hsv(nAng);
for l0 = 1:nAng
plot(TPfin(:,:,l0).','Color',c(l0,:),'LineWidth',1)
end
plot(rigidDModes,'ko','MarkerFaceColor','k')%,'MarkerSize',8)
plot(rigidAModes,'ks','MarkerFaceColor','k')%,'MarkerSize',8)
% Convected mode
cMode = -omega*newADData.delta;
plot(real(cMode),imag(cMode),'k^','MarkerFaceColor','k')%,'MarkerSize',25)
%set(gca,'color',.9*[1,1,1]);
axis equal

% plot(smallRootsP,'g*')
% plot(smallRootsM,'b*')


hold off
axis([-15,10,-17,17])
xlabel('$\Re[\theta]$');%,'Interpreter','Latex')
ylabel('$\Im[\theta]$');%,'Interpreter','LaTeX')
set(gca,'Layer','top')

% cleanfigure
% matlab2tikz([imageFolder,'mode-traj-',num2str(m),'.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
%                  'extraaxisoptions','ylabel shift = -.4cm,,axis on top');
%%             
% if m==2 
figure(2)


Csmall = logspace(-4,0);
Clarge = logspace(-1,3);
nAsym = (0:(nDM-1)).';
dN = ones(nDM,1); dN(1)=2;
zeta = @(zVar) mysqrt(omega*w,zVar);
smallRootsP = rDP + Csmall.*1i.*(rDP-cMode)./(dN.*sPG.*rDP).*(1-(-1).^nAsym.*cos(dPG.*rDP+sigma));
smallRootsM = rDM + Csmall.*1i.*(rDM-cMode)./(dN.*sPG.*rDM).*(1-(-1).^nAsym.*cos(dPG.*rDM+sigma));

zP = (dPG*rAP + delPG*fm)/sPG;
zM =-(dPG*rAM + delPG*fm)/sPG;

largeRootsP = rAP + 1./Clarge.'.*1i.*(zP.^2)./((rAP-cMode).*delPG.*zeta(fm));
largeRootsM = rAM - 1./Clarge.'.*1i.*(zM.^2)./((rAM-cMode).*delPG.*zeta(fm));
largeCmode = cMode+1./Clarge*1i*zeta(cMode).*sin(sPG*zeta(cMode))./(cos(sPG*zeta(cMode))-cos(dPG*cMode+sigma));
recCol1 = [255,240,240]/255;
recCol2 = [240,240,255]/255;

semilogx(abs(C1vec(:,1)),imag(TPfin(:,:,1)),'k','LineWidth',2)
hold on
semilogx(Csmall,imag(smallRootsP).','b--')
semilogx(Clarge,imag(largeCmode).','r--')
semilogx(Clarge,imag(largeRootsP(:,[1,4,5,6,7,8])),'r--')
hold off
axis([1e-4,1e2,-.01,18])
xlabel('$C_{II}$');%,'Interpreter','Latex')
ylabel('$\Im[\theta]$');%,'Interpreter','LaTeX')
set(gca,'color',recCol1);

% cleanfigure
% matlab2tikz([imageFolder,'asympU.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
%                  'extraaxisoptions','axis on top');


figure(3)
semilogx(abs(C1vec(:,1)),imag(TPfin(:,:,1)),'k','LineWidth',2)
hold on
semilogx(Csmall,imag(smallRootsM).','b--')
semilogx(Clarge,imag(largeCmode).','r-.')
semilogx(Clarge,imag(largeRootsM(:,[2,5,6,7,8,11])),'r--')
hold off
axis([1e-4,1e2,-18,0.01])
xlabel('$C_{II}$');%,'Interpreter','Latex')
ylabel('$\Im[\theta]$');%,'Interpreter','LaTeX')
set(gca,'color',recCol2);
1
% cleanfigure
% matlab2tikz([imageFolder,'asympL.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
%                  'extraaxisoptions','axis on top');

end
