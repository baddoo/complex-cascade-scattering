%% Makes a matrix periodic

function M = funPer(mat,change,ns)

    n = -ns:ns;
    n3 = permute(n,[1,3,2]);

    matSize = [size(mat,1),size(mat,2)*(2*ns+1)];
    M=reshape(mat+n3*change,matSize);

end