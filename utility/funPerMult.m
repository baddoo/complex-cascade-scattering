function M = funPerMult(mat,change,ns)
%FUNPERMULT Makes a matrix periodic but with a multiplicative change between period
% windows
n = -ns:ns;
n3 = permute(n,[1,3,2]);

matSize = [size(mat,1),size(mat,2)*(2*ns+1)];
M=reshape(mat.*(change.^n3),matSize);

end