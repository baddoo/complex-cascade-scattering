% Tests whether the asymptotic approximations for the roots are correct.
AAData.w = 3;
d = .9; s = .5;
ADData.spac = [s,d];
ADData.chie = atan(d/s);
AAData.Sigma = 1; AAData.omega = 10;
Modes.trunc = 100;

tol = 1e-3;
for j = 1:3

    figure(j)
    clf
    hold on
switch j
    case 1
        ADData.mu = [rand + 1i*rand 0 0];  
        title('asymptotic approximation errors for case I boundary condition','interpreter','latex','fontsize',20)
    case 2
        ADData.mu = [rand+1i*rand rand+1i*rand 0];
        title('asymptotic approximation errors for case II boundary condition','interpreter','latex','fontsize',20)
    case 3
       ADData.mu = [rand+1i*rand rand+1i*rand rand+1i*rand];   
       title('asymptotic approximation errors for case III boundary condition','interpreter','latex','fontsize',20)

end
    
f = @(xVar) 1./KNumlogD(xVar,ADData,AAData);
logD = @(xVar) KNumlogD(xVar,ADData,AAData);

ADData.case = j;
asympGuess = computeAsympGuess(ADData,AAData,Modes);
knownRootsInit = newtonKRoots(asympGuess,f,logD,tol,[]);

% Divide into quadrants
err = abs(asympGuess - knownRootsInit);
err1 = err(0<angle(asympGuess) & angle(asympGuess)<pi/2);
err2 = err(pi/2<angle(asympGuess) & angle(asympGuess)<pi);
err3 = err(0>angle(asympGuess) & angle(asympGuess)>-pi/2);
err4 = err(-pi<angle(asympGuess) & angle(asympGuess)<-pi/2);

q1 = semilogy(err1,'-','LineWidth',3);
q2 = semilogy(err2,'-','LineWidth',3);
q3 = semilogy(err3,'-','LineWidth',3);
q4 = semilogy(err4,'-','LineWidth',3);

grid on

ylabel('error','Interpreter','latex','FontSize',20)
xlabel('$n$','Interpreter','latex','FontSize',20)

set(gca,'YScale','log')

legend('quadrant 1','quadrant 2','quadrant 3','quadrant 4')


hold off


end
