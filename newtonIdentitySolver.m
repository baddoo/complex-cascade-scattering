% Finds the roors of the polynomials satisfying Newton's identities
%
% See the wikipedia page on Newton's identities for further details.
% p are the Newton sums
% order is the order of the highest-order sum
% nCases is the number of cases -- i.e. how many annuli are considered.
% maxOrder is a vector that denotes the maximum order of each case.

function rootVec = newtonIdentitySolver(p,maxOrder)

order = size(p,1);
nCases = size(p,2);
e = ones(order+1,nCases);

for l = 1:order
     j = (1:l).';
     size(e);
     sum((-1).^(j-1).*e(l-j+1,:).*p(1:l,:),1);
     e(l+1,:) = 1./l.*sum((-1).^(j-1).*e(l-j+1,:).*p(1:l,:),1);
end

for i = 1:nCases
e((maxOrder(i))+2:end,i) = 0;
end

rootMat = (-1).^(1:order+1)'.*e;
rootVec = cell2mat(cellfun(@roots, num2cell(rootMat, 1), 'UniformOutput', false));

end