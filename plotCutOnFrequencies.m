function [RmodesCO,AmodesCO]=plotCutOnFrequencies(Cp,newADData,newAAData,Modes)
% Calculate the cut on frequencies of various modes
sigmao=newAAData.sigmao;h=newADData.spac(1); d=newADData.spac(2);
s=sqrt(h^2+d^2); k=newAAData.k;
delta=newADData.delta; M=newADData.M;
w=newAAData.w;

Rm=((-Modes.amodes):(Modes.amodes)).';
RmodesCO=[(sigmao-2*Rm*pi)/(w*s-d*delta*M^2);...
          -(sigmao-2*Rm*pi)/(w*s+d*delta*M^2)];
RmodesCO(RmodesCO<min(k) | RmodesCO>max(k))=[];
[~,indR]=min(abs(bsxfun(@minus,k.',RmodesCO.')));


Am=(-Modes.dmodes:Modes.dmodes).';
AmodesCO=-Am*pi/(w*h);

[~,indA]=min(abs(bsxfun(@minus,k.',AmodesCO.')));

m=-1; % Doesn't seem to do anything
GModesMeetDModes=[-pi/(h*(w.^2-delta.^2))*(2*h*m*delta)+sqrt((2*h*m*w).^2+Am.^2*(w^2-delta.^2));...
                  -pi/(h*(w.^2-delta.^2))*(2*h*m*delta)-sqrt((2*h*m*w).^2+Am.^2*(w^2-delta.^2))];
[~,indAGD]=min(abs(bsxfun(@minus,k.',GModesMeetDModes.')));

hold on
ColourSpecification;
plot(k(indR)/2,abs(Cp(indR)),'x','MarkerFaceColor',co(2,:),'MarkerEdgeColor',co(2,:)) %#ok<IDISVAR,NODEF>
plot(k(indA)/2,abs(Cp(indA)),'x','MarkerFaceColor',co(3,:),'MarkerEdgeColor',co(3,:))
%plot(GModesMeetDModes,abs(Cp(indAGD)),'o','MarkerFaceColor',co(4,:),'MarkerEdgeColor','k')
hold off
end