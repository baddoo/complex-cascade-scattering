function [TP3,TM3,asympGuess] = findDuctModes2(R0,chi,ADData,AAData,Modes)

tol = 1e-10;

f = @(xVar) 1./KNumlogD(xVar,ADData,AAData);
logD = @(xVar) KNumlogD(xVar,ADData,AAData);

%if mu == zeros(size(ADData.mu))
asympGuess = computeAsympGuess(ADData,AAData,Modes);
knownRootsInit = newtonKRoots(asympGuess,f,logD,tol,[]);
knownRootsInit = myUnique(knownRootsInit,tol);

% Only use the roots that are relatively close to the domain
loc = abs(knownRootsInit)<2*R0;
allRootsUnique = rootFinder(f,logD,knownRootsInit(loc),R0,chi,tol);
allRootsUnique = [allRootsUnique;knownRootsInit(~loc).'];

% Separate into components in the upper and lower half planes
TP = sort(allRootsUnique(imag(allRootsUnique)>=0));
TM = sort(allRootsUnique(imag(allRootsUnique)<0));

% Reshape into third dimension
TP3 = permute(TP,[3,2,1]);
TM3 = permute(TM,[3,2,1]);

end