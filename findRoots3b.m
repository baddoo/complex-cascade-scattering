function F = findRoots3b(s,d,sigma,mu,k,w,x0)

options=optimset('Display','off');%,'LargeScale','off','TolFun',1e-8,'MaxIter',1e5,'MaxFunEvals',1e5);
%options = [];
[F,~,exitflag] = fsolve(@rootFun,x0,options);
if exitflag < 1
  %F = [];
end

%     function G =rootFunB(x)
%         G= (mysqrt(k*w,x).*sin(s*mysqrt(k*w,x)) + (mu)*(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma))) ...
%             +0./(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma));
%     end

    function G =rootFunB(x)
          G= mysqrt(k*w,x).*sin(s*mysqrt(k*w,x))./(x.*(cos(s*mysqrt(k*w,x)) - cos(d*x + sigma))) ...
              + mu(3)./x;% +mu(2).*x + mu(1);
    end

    function H = rootFun(x)
    H(1,:) = real(rootFunB(x(1,:)+1i*x(2,:)));
    H(2,:) = imag(rootFunB(x(1,:)+1i*x(2,:)));
    end

end