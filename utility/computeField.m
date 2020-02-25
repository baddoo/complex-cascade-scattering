function acousticField = computeField(data,type)

X=data.X;
Y=data.Y;
s = data.spac(1); d = data.spac(2); kx = data.kx;


upstream = u(data,type);
upperTriangle = acousticFieldUpperTriangle(data,type);
rectangle = acousticFieldRectangle(data,type);
lowerTriangle = acousticFieldLowerTriangle(data,type);
downstream = acousticFieldDownstream(data,type);

m = [1,1,1,1,size(kx,5)];

% Define complexified coordinate
Z = X + 1i*Y;

h = @(zVar) evaluateField(zVar,data);

end

% upstream(repmat(X>=Y*d/s,m))=0;
% upperTriangle(repmat(X<Y*d/s | X>=d,m))=0;
% rectangle(repmat(X<d | X>=2,m))=0;
% lowerTriangle(repmat(X>=2+Y*d/s | X<2,m))=0;
% downstream(repmat(X-2<Y*d/s,m))=0;
%

