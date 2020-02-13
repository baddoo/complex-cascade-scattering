% custom newton root finder

function [x,steps] = myNewton(f,fp,x0,tol)

err = inf;
xn = x0;
j=0;

while err>tol
    fn=f(xn); %Calculating the value of function at x0
    fn_der=fp(xn); %Calculating the value of function derivative at x0
    xn=xn-fn./fn_der; % The Formula
    err = norm(f(xn),'inf');
    j = j+1;
end

steps = j;
x = xn;

end