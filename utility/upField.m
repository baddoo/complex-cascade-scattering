function phi = upField(Z,data,type)

X = real(Z); Y = imag(Z);

del=data.spac(3);
SQRTa=data.SQRTa;
LPa = data.LPa;
ZPa = data.ZPa;
data.comb=[1,0,1,0];

%Ufin= pi*1i*sum(bsxfun(@times,permute(D(LPa,data),[3,2,1,4,5]),(A1aResP+A1bResP)),3);
phi = pi*1i*sum(bsxfun(@times,permute(D(LPa,data),[3,2,1,4,5]),(A1aResP+A1bResP)),3);

%Dcoefs = permute(D(permute(LPa,[3,2,1]),data),[3,2,1]);
%phi = pi/del*sum(Dcoefs.*ZPa./SQRTa.*exp(-1i*(X.*LPa - Y.*ZPa)),3);

end