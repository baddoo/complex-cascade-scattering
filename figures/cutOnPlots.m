imageFolder = '../../images/';


del = .8;
M = .4;
chi = [10,30,30]*pi/180;

nPts = 2000;
nDM = 2;
nAng = 4;
TPfin = zeros(2*nDM,nPts,nAng);

xSol = zeros(2*nDM,1);
C1vec = [0,0.1,.5,1];%exp(1i*ang);
omegaVec = linspace(0.001,20,nPts);
    
for l0 = 1:nAng
    for l1 = 1:nPts
        C1 = C1vec(l0);
ADData=struct('spac',   [del*cos(chi(1)),del*sin(chi(1))],...%Blade Spacing [h,d]
              'physC',   1,...                          %Blade length
              'c0',      340,...                        %Speed of sound
              'M',       M,...                         %Mach number 
              'case',    'II', ...                       %Case
              'C',       C1);                            %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omega',    omegaVec(l1)*(1+1e-9i),...                      %Tangential Frequency
               'kn',       [],...                       %Normal frequency
               'k',        2*(1+1e-9i),...
               'kz',       0,...                       %Spanwise frequency
               'sigmao',   4*pi/3*(1+1e-9i),...                  %Interblade phase angle in (x,y)-space
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
  
   options=optimset('Display','off','LargeScale','off','TolFun',1e-16,'MaxIter',10000,'MaxFunEvals',10000);
    for j = 1:(2*nDM)
    xSol(j) = fsolve(F,TPfin(j,l1-1,l0));  
    end
    TPfin(:,l1,l0) = xSol;
    
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
rectangle('Position',[-6.5,0,13,6.5],'FaceColor',recCol1,'EdgeColor',recCol1)
rectangle('Position',[-6.5,-6.5,13,6.5],'FaceColor',recCol2,'EdgeColor',recCol2)
% probLoc = find(ang==pi*3/2);
% loc1 = find(abs(TPfin(1,:,probLoc))>10);
% %loc2 = find(abs(TPfin(2,:,probLoc))>10);
% if ~isempty(loc1); TPfin(1,loc1(1):end,probLoc)=nan; end
%TPfin(2,loc2(2):end,probLoc)=nan;

c = hot(ceil(1.5*nAng));
for l0 = 1:nAng
plot(TPfin([2,4],2:end,l0).','-','Color',c(l0,:),'LineWidth',1)
plot(TPfin([2,4],2,l0).','o','Color',c(l0,:),'MarkerFaceColor',c(l0,:))

end
%plot(rigidDModes,'ko','MarkerFaceColor','k')%,'MarkerSize',8)
%plot(rigidAModes,'ks','MarkerFaceColor','k')%,'MarkerSize',8)
% Convected mode
%cMode = -omega*newADData.delta;
%plot(real(cMode),imag(cMode),'k^','MarkerFaceColor','k')%,'MarkerSize',25)
%set(gca,'color',.9*[1,1,1]);
axis equal
hold off
% if m==1
% axis([-10,10,-20,20])
% elseif m==2
%     axis([-12,12,-20,20])
% elseif m==3
%     axis([-12,12,-20,20])
% end
axis([-6,6,-6,6])
xlabel('$\Re[\theta]$');%,'Interpreter','Latex')
ylabel('$\Im[\theta]$');%,'Interpreter','LaTeX')
box on
%set(gca,'Layer','top')

cleanfigure
matlab2tikz([imageFolder,'mode-traj-freq.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
                  'extraaxisoptions','ylabel shift = -.2cm,axis on top');
%%             
% if m==2             
figure(2)
clf
for l0 = 1:nAng
plot(abs(omegaVec(2:end)),abs(imag(TPfin(2,2:end,l0))).','Color',c(l0,:),'LineWidth',2)
hold on

end
set(gca,'color',recCol1);
ylabel('$\Im[\theta]$');%,'Interpreter','Latex')
xlabel('$\omega$');%,'Interpreter','LaTeX')
hold off
cleanfigure
 matlab2tikz([imageFolder,'im-part-1.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
                  'extraaxisoptions','ylabel shift = -.1cm');

figure(3)
clf
for l0 = 1:nAng
plot(abs(omegaVec(2:end)),(imag(TPfin(4,2:end,l0))).','Color',c(l0,:),'LineWidth',2)
hold on

end
set(gca,'color',recCol2);
ylabel('$\Im[\theta]$');%,'Interpreter','Latex')
xlabel('$\omega$');%,'Interpreter','LaTeX')
hold off
ylim([-inf,eps+.01])
cleanfigure
 matlab2tikz([imageFolder,'im-part-2.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
                  'extraaxisoptions','ylabel shift = -.1cm');

      