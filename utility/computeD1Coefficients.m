function D1Data= computeD1Coefficients(coefdata)

extractStructureData;
coefdata.AAData;
w0  = An;
T  = w0/(4i*pi^2*KMGM0);

S = 0;

V = w0*(PM0-GM0)./(4*pi^2*KPGM0);

Tsum = permute(sum(T./TMminGM0,3),[1,2,4,3,5]);

A = (1i*(TMd-PM0)./KpprTM).*Tsum;
A1 = permute(A,[3,4,5,1,2]);
V1 = permute(V,[3,4,5,1,2]);

B1=zeros(size(A1));
C1=zeros(size(A1));
Dres = zeros(size(A1));

% Can remove this looping since it is now obselete
for ifreq = 1:nfreq
Dres(:,:,ifreq) = F1(:,:,ifreq)*A1(:,:,ifreq) + Vpreq(:,:,ifreq)*V1;
B1(:,:,ifreq)=invmat(:,:,ifreq)*Dres(:,:,ifreq);
C1(:,:,ifreq)=L1(:,:,ifreq)*B1(:,:,ifreq);
end

B=permute(B1,[5,4,1,2,3]);
C=permute(C1,[5,4,1,2,3]);

D1Data.T =  T;
D1Data.S =  S;
D1Data.A  = A;
D1Data.V  = V;
D1Data.B  = B;
D1Data.C  = C;

end
