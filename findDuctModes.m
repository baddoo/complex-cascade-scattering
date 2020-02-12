function [TP3,TM3,asympGuess] = findDuctModes(ADData,AAData,Modes)


[asympRoots,asympGuess,boxHeight] = findAsympRoots(ADData,AAData,Modes);

closeRoots = findCloseRoots(boxHeight,ADData,AAData);
allRoots = [closeRoots,asympRoots];
tol = 1e-10;
% Remove non-unique entries
allRootsUnique = myUnique(allRoots,tol);

% Separate into components in the upper and lower half planes
TP = sort(allRootsUnique(imag(allRootsUnique)>0));
TM = sort(allRootsUnique(imag(allRootsUnique)<0));

% Reshape into third dimension
TP3 = permute(TP,[1,3,2]);
TM3 = permute(TM,[1,3,2]);

end