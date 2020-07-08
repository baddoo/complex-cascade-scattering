function x = newtonKRoots(x0,f,logD,tol,knownRoots)
%newtonKRoots  Finds the roots of a meromorphic function using Newton
%iteration.
%X = newtonKRoots(X0,F,LOGD,TOL,KNOWNROOTS) returns the roots of the function F
% found by Newton iteration at initial point X0 to a tolerance of TOL. The
% logarithmic derivative (F'/F) is supplied as LOGD. Any previously known
% roots may also be supplied as KNOWNROOTS to ensure that the algorithm
% converges to new roots.

err = inf;
xn = x0;
j=0;

% Define new logarithmic derivative with known roots removed. Uses the 5-th
% dimension for historical reasons.
polyLogD = @(xVar,krV) sum(1./(xVar - permute(knownRoots(:),[5,2,3,4,1])),5);

% Remove nan and entries with infinite gradient
locN= isnan(logD(x0)) | logical(logD(x0)==0);
xn = xn(~locN);

while err>tol 
    xn=xn-1./(logD(xn)-polyLogD(xn)); % The Newton iteration formula
    xn(isnan(f(xn)))=[]; % Remove entries where the function becomes undefined.
    err = norm(f(xn),'inf'); % Calculate error
    j = j+1;
end

x = xn;

end