function acousticField = computeField(data,type)

X=data.X;
Y=data.Y;
he = data.spac(1); d = data.spac(2); k5 = data.k5;

upstream = acousticFieldUpstream(data,type);
upperTriangle = acousticFieldUpperTriangle(data,type);
rectangle = acousticFieldRectangle(data,type);
lowerTriangle = acousticFieldLowerTriangle(data,type);
downstream = acousticFieldDownstream(data,type);

m = [1,1,1,1,size(k5,5)];

upstream(repmat(X>=Y*d/he,m))=0;
upperTriangle(repmat(X<Y*d/he | X>=d,m))=0;
rectangle(repmat(X<d | X>=1,m))=0;
lowerTriangle(repmat(X>=1+Y*d/he | X<1,m))=0;
downstream(repmat(X-1<Y*d/he,m))=0;

acousticField =  upstream + upperTriangle + rectangle + lowerTriangle + downstream; 

end