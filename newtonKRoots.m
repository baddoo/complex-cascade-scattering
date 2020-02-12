% Find the roots of the Wiener--Hopf kernel, K, with regularised functions.

function x = newtonKRoots(x0,f,logD,tol,knownRoots)

err = inf;
xn = x0;
j=0;

polyLogD = @(xVar,krV) sum(1./(xVar - permute(knownRoots(:),[5,2,3,4,1])),5);

% Remove nan and entries with infinite gradient
locN= isnan(logD(x0)) | logical(logD(x0)==0);
xn = xn(~locN);

while err>tol 
    xn=xn-1./(logD(xn)-polyLogD(xn)); % The Formula
    xn(isnan(f(xn)))=[]; % Remove entries where the function becomes undefined.
    err = norm(f(xn),'inf'); % Calculate error
    j = j+1;
end

x = xn;

end