function h = evaluateField(Z,data,type)

% Split up Z
X = real(Z); Y = imag(Z);
s=data.spac(1); d=data.spac(2);

% Partition Z into relevant domains
upInd = find(X<=Y*d/s);
upIntInd = find(X>Y*d/s & X<=d); % Z in the upper inter-blade region
intInd = find(X>d & X<=2); % Z in the inter-blade region
dwnIntInd = find(X<=2+Y*d/s & X>2); % Z in the lower inter-blade region
dwnInd = find(X-2>Y*d/s); % Z in the downstream

h = zeros(size(Z));
h(upInd) = upField(Z(upInd),data,type);
h(upIntInd) = upIntField(Z(upIntInd),data,type);
h(intInd) = intField(Z(intInd),data,type);
h(dwnIntInd) = dwnIntField(Z(dwnIntInd),data,type);
h(dwnInd) = upField(Z(dwnInd),data,type);

%  h = reshape([upField(upZ,data,type); ...
%               upIntField(upIntZ,data,type); ...
%               intField(intZ,data,type); ...
%               dwnIntField(dwnIntZ,data,type); ...
%               dwnField(dwnZ,data,type)],size(Z));

end