imageFolder = '../../../LaTeX/porous-jfm-r1/images/';

chordDim = 1;
vaneDistDim = .8*chordDim;
M = .4;
chi = 10*pi/180;

nPts = 2000;
nDM = 2;
nAng = 4;
TPfin = zeros(2*nDM,nPts,nAng);

xSol = zeros(2*nDM,1);
C = [1e-5,0.1,.5,1];
omegaVec = linspace(0.001,20,nPts);
    
for l0 = 1:nAng
    for l1 = 1:nPts
        C1 = C(l0);
ADData=struct('spacDim',   [vaneDistDim*cos(chi),vaneDistDim*sin(chi)],...%Blade Spacing [h,d]
              'chordDim',   1,...                          %Blade length
              'c0',         340,...                        %Speed of sound
              'M',          M,...                         %Mach number 
              'Wdim',       0,...
              'case',       2, ...                       %Case
              'C',          C1);                            %Coefficient of case                       
             
%% Aeroacoustic Data
AAData=struct( 'omegaDim', [],...
               'omega',    omegaVec(l1),...                      %Tangential Frequency
               'kx',       2,...
               'kxDim',    [],...
               'ky',       [],...
               'kyDim',    [],...
               'kz',       0,...                       %Spanwise frequency
               'kzDim',    [],...                       %Spanwise frequency
               'sigmao',   4*pi/3,...                  %Interblade phase angle in (x,y)-space
               'Amp',      [nan,1,nan]); %Amplitude of gust in form [At,An,A3]
% Prepare data
[newADData,newAAData] = prepareData(ADData,AAData);    

    mu = newADData.mu;
    dPG = newADData.spac(2); sPG = newADData.spac(1); omega = newAAData.omega;
    delPG = newADData.spac(3);
    w = newAAData.w;
    sigma = newAAData.sigma;
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
    tol = 1e-8;
    TPfin(:,l1,l0) = newtonKRoots(TPfin(:,l1-1,l0),f,logD,tol,[]);   
    
    end
[l0,l1]
    end
end


%% Plots
figure(1)

% Remove misleading point
loc = find(imag(TPfin([2],:,1))>1e-3,1,'last');
TPfin([2,4],loc+1,1)=0;
recCol1 = [255,240,240]/255;
recCol2 = [240,240,255]/255;
rectangle('Position',[-6.5,0,13,6.5],'FaceColor',recCol1,'EdgeColor',recCol1)
hold on
rectangle('Position',[-6.5,-6.5,13,6.5],'FaceColor',recCol2,'EdgeColor',recCol2)

c = hot(ceil(1.3*nAng));
plot(TPfin([2],2:end,1).','-','Color','k','LineWidth',3)
legendInfo{1} = ['$C_{II} = ' num2str(0) '$']; 
%plot(TPfin([2,4],2,l0).','o','Color',c(l0,:),'MarkerFaceColor',c(l0,:))

for l0 = 1:(nAng-1)
    plot(TPfin([2],2:end,l0+1).','-','Color',c(l0+1,:),'LineWidth',2)
    legendInfo{l0+1} = ['$C_{II} = ' num2str(C(l0+1)) '$']; 

%length = sum(abs(diff(TPfin([2,4],:,l0),1,2))); % Calculate length of each trajectory
% %[~,loc] = min(abs(cumsum(abs(diff(TPfin(j,:,l0),1,2)))-length/2));
% if l0 == 1
% loc2 = 21;
% arPts2 = [TPfin(2,loc2,l0),TPfin(2,loc2+1,l0)]; arPts2= [real(arPts2);imag(arPts2)];
% arrow(arPts2(:,1),arPts2(:,2),10,'FaceColor',c(l0,:),'EdgeColor',c(l0,:))
% loc4 = 21;
% arPts4 = [TPfin(4,loc4,l0),TPfin(4,loc4+1,l0)]; arPts4= [real(arPts4);imag(arPts4)];
% arrow(arPts4(:,1),arPts4(:,2),10,'FaceColor',c(l0,:),'EdgeColor',c(l0,:))
% end
% 
% if l0 ==3
% loc2 = 35;
% arPts2 = [TPfin(2,loc2,l0),TPfin(2,loc2+1,l0)]; arPts2= [real(arPts2);imag(arPts2)];
% arrow(arPts2(:,1),arPts2(:,2),10,'FaceColor',c(l0,:),'EdgeColor',c(l0,:))
% loc4 = 20;
% arPts4 = [TPfin(4,loc4,l0),TPfin(4,loc4+1,l0)]; arPts4= [real(arPts4);imag(arPts4)];
% arrow(arPts4(:,1),arPts4(:,2),10,'FaceColor',c(l0,:),'EdgeColor',c(l0,:))
% end

end

plot(TPfin([4],2:end,1).','-','Color','k','LineWidth',3)
plot(TPfin([2,4],2,1).','o','Color','k','MarkerFaceColor','k')

for l0 = 1:(nAng-1)
    plot(TPfin([4],2:end,l0+1).','-','Color',c(l0+1,:),'LineWidth',2)
    plot(TPfin([2,4],2,l0+1).','o','Color',c(l0+1,:),'MarkerFaceColor',c(l0+1,:))
end
legend(legendInfo,'Interpreter','latex','Location','northwest')
hold off

axis equal
hold off

axis([-6,6,-3,3])
xlabel('$\Re[\theta]$');
ylabel('$\Im[\theta]$');
box on

cleanfigure
matlab2tikz([imageFolder,'mode-traj-freq.tex'], 'height', '\fheight', 'width', '\fwidth','parseStrings',false,'extratikzpictureoptions','trim axis left, trim axis right',...
                  'extraaxisoptions','ylabel shift = -.2cm,axis on top');
%%             
% if m==2             
figure(2)
clf
plot(abs(omegaVec(2:end)),abs(imag(TPfin(2,2:end,1))).','Color','k','LineWidth',3)
hold on
for l0 = 1:(nAng-1)
plot(abs(omegaVec(2:end)),abs(imag(TPfin(2,2:end,l0+1))).','Color',c(l0+1,:),'LineWidth',2)
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
plot(abs(omegaVec(2:end)),(imag(TPfin(4,2:end,1))).','Color','k','LineWidth',2)
hold on
for l0 = 1:(nAng-1)
plot(abs(omegaVec(2:end)),(imag(TPfin(4,2:end,l0+1))).','Color',c(l0+1,:),'LineWidth',2)
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

      