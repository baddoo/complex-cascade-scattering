function D3=D3(data,data2)

B=data2.B;
gminTP=data.gminTP;
KPG=data.KPG;
KPTP=data.KPTP;

num=B.*KPTP;
denom=gminTP;

D3terms=bsxfun(@rdivide,num,denom);

D3=-sum(D3terms,3)./KPG;

end