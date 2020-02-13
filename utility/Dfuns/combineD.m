function D=combineD(data,dnum)

if data.comb(1)==1; D_O1 = D1(data,dnum); else D_O1=0; end %#ok<*SEPEX>
if data.comb(2)==1; D_O2 = D2(data,dnum); else D_O2=0; end
if data.comb(3)==1; D_O3 = D3(data,dnum); else D_O3=0; end
if data.comb(4)==1; D_O4 = D4(data,dnum); else D_O4=0; end

D=D_O1+D_O2+D_O3+D_O4;

end