%% Removes duplicate entries in a complex vector

function finV = myUnique(v,tol)

vr = real(v(:));
vi = imag(v(:));

C = uniquetol([vr,vi],tol,'ByRows',true,'DataScale',max(abs(v(:))));

finV = reshape(C(:,1)+1i*C(:,2),size(v(1:numel(C)/2)));

end