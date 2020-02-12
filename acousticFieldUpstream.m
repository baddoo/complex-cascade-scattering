function acousticFieldUpstream = acousticFieldUpstream(data,type)

A1aResP = data.A1aResP.(type);
A1bResP = data.A1bResP.(type);
LPa = permute(data.LPa,[3,2,1,4,5]);

data.comb=[1,0,1,0];

Ufin= pi*1i*sum(bsxfun(@times,permute(D(LPa,data),[3,2,1,4,5]),(A1aResP+A1bResP)),3);

acousticFieldUpstream = Ufin;

end