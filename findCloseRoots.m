function closeRoots = findCloseRoots(boxHeight,ADData,AAData)

mu = ADData.mu;
he = ADData.spac(1);
d = ADData.spac(2);
sigma5 = AAData.sigma5;
w = AAData.w;
omega5 = AAData.omega5;
disp(['The box height is ' num2str(round(boxHeight))]);
%boxHeight= 1*boxHeight
zeta = @(xVar) mysqrt(omega5*w,xVar);

fn  = @(x) zeta(x).*sin(he*zeta(x))+ (mu(1)-1i*x*mu(2)-x.^2*mu(3)).*(cos(he*zeta(x)) - cos(d*x + sigma5));

fnp = @(x) -x./zeta(x).*sin(he*zeta(x))-he*x.*cos(he*zeta(x)) ...
           + (mu(1)-1i*x*mu(2)-x.^2*mu(3)).*(he*x./zeta(x).*sin(he*zeta(x)) + d*sin(d*x+sigma5)) ...
           + (-1i*mu(2)-2*x*mu(3)).*(cos(he*zeta(x)) - cos(d*x + sigma5));

argFun = @(xVar) fnp(xVar)./fn(xVar);

boxLength=60;

ways = [(boxLength-1i*boxHeight)/2,(boxLength+1i*boxHeight)/2,(-boxLength+1i*boxHeight)/2];
Npr = (real(integral(argFun,-(boxLength+1i*boxHeight)/2,-(boxLength+1i*boxHeight)/2,'Waypoints',ways,'ArrayValued',true)/(2i*pi)));

N = round(Npr);
N3 = permute(1:N,[1,3,2]);

polyFun = @(xVar) bsxfun(@power,xVar,N3);
newFun = @(zVar) polyFun(zVar).*argFun(zVar);

p = integral(newFun,-(boxLength+1i*boxHeight)/2,-(boxLength+1i*boxHeight)/2,'Waypoints',ways,'ArrayValued',true)/(2i*pi);

e = ones(1,N+1);

for l = 1 :N
    j = (1:l).';
    e(l+1) = 1./l.*sum(permute((-1).^(j-1).'.*e(l-j+1),[1,3,2]).*p(:,:,(1:l)),3);
end

closeRootsPrelim = roots((-1).^(1:N+1).*e);
closeRoots = zeros(size(closeRootsPrelim.'));

for nL = 1:numel(closeRoots)
 closeRoots(nL) = findRoots3(he,d,sigma5,mu,omega5,w,closeRootsPrelim(nL));
end

end

