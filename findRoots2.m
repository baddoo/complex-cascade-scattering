function xSol = findRoots2(s,d,sigma,mu,k,w,x0)

options=optimset('Display','off','LargeScale','off','TolFun',1e-16,'MaxIter',2000,'MaxFunEvals',10000);
[xSol,~,exitflag] = fsolve(@rootFunB,x0,options);
if exitflag < 1
 % xSol = [];
end
    function F =rootFunB(x)
        %F= mysqrt(k*w,x).*sin(s*mysqrt(k*w,x)) + mu*(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma))...
        %    ./(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma));
          F= mysqrt(k*w,x).*sin(s*mysqrt(k*w,x))./(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma)) + mu(1);
    end

%     function F = rootFun(x)
%     F(1,:) = real(rootFunB(x(1,:)+1i*x(2,:)));
%     F(2,:) = imag(rootFunB(x(1,:)+1i*x(2,:)));
%     end

end