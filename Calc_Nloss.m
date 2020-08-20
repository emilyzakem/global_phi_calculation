
%add up the volume of water with phi at 1 or below:
anvolall=[];
Nlossall=[];
threshhall=0:0.1:2;
%for j=1:length(threshhall);
%threshh=threshhall(j); %MOST sensitive to this!
threshh=1;
%thresho=5; %5 uM O2



phib=phi;
phib(phib>=threshh)=nan;
phib(~isnan(phib))=1;

%phib(o2conc>10)=nan;



o2b=o2conc;
o2b(o2b>threshh)=nan;
o2b(~isnan(o2b))=1;

dzfgrid = repmat(dzf,1,size(poc3di,2),size(poc3di,3));
dxgrid=repmat(dx,1,size(phi,1),size(phi,3));
dxgrid=permute(dxgrid,[2 1 3]);

vol = dzfgrid.*dy.*dxgrid; %m3


%anaerobic volume in km3
anvol=nansum(nansum(nansum(phib.*vol)))*(1e-3)^3*1e-6 %from m3 to 1e6 km3 of water
anvolall(j)=anvol;
o2vol=nansum(nansum(nansum(o2b.*vol)))*(1e-3)^3*1e-6 %from m3 to 1e6 km3 of water


%how much denitrification (N loss) should be happening here, then?
%in each box where phi<1, take the diff(POC flux) to get moles C /day
%POC in: mmol C/m3/day (molC/m2/day/dzc)

%Three ways of calculating this for one phi, which is set above 
%(phib is either from phi or from phidiv): which of the 3 is best?
%1. 235 for phi (30 for phidiv)
%POCb=POCin.*phib.*vol; %multiplied by phi mask and volume: mmoles C/day
%POCbdiff=POCb(2:end,:,:)-POCb(1:end-1,:,:); %left in each box
%2. 415
%this way gives more- almost double (bc i don't cut off the 2:end edges):
%POCbdiff=(POCin(2:end,:,:)-POCin(1:end-1,:,:)).*phib(1:21,:,:).*vol(1:21,:,:); %left in each box
%3. 302
%could use POCdiv, and do the following: it gives a different number: 302 instead of 235.
%i need to check the way depth is applied.. i THINK this is right bc the
%depth of just the one box is applied to the divergence. only coincidentally does it
%also give the middle value here.
POCbdiff=POCdiv.*phib.*vol; %left in each box

%now get yield for denitr: look at choices of Bianchi 2012:
%but go back to Gruber and Sarmiento for still the 'best' estimate:
yod=.2;
ynd=yod*.8; %always 20% of aerobic
    yne=ynd/20*5/(1-ynd)*5;
rN=ynd/yne;%my estimate is low!

rN=104/106; %Gruber and Sarmiento 1997: 104+-15 : 106
%rN = (104+15)/106;

Nloss3D=POCbdiff*rN; %mmoles NO3 per day needed to oxidize all that POC
Nloss2D = squeeze(nansum(Nloss3D,1));
Nloss2D(Nloss2D == 0) = nan;
Nloss = nansum(nansum(nansum(Nloss3D))); %N lost to C oxidation
Nloss = Nloss*14*1e-3*1e-12*365 %Tg N/year (from mmoles N *14g/mol * 1e-3moles/mmol * 1e-12 Tg/g * 365 days/year
%Nlossall(j)=Nloss;


%%
%Feb 2017:
%try to calculate "residual" OM in this domain, using the info of how much
%O2 has been consumed:
%POCbdiff = 
%O2bdiff = posO2div.*phib.*vol; %this is the total O2 within phi. from O2div within phi
%O2bdiff = O2div.*phib.*vol; %this is the total O2 within phi. from O2div within phi

%residPOC = POCbdiff - O2bdiff/r; %how much "leftover" POC?
%Nloss=residPOC*rN; %mmoles NO3 per day needed to oxidize all that POC
%Nloss = nansum(nansum(nansum(Nloss))); %N lost to C oxidation
%Nloss = Nloss*14*1e-3*1e-12*365 %Tg N/year (from mmoles N *14g/mol * 1e-3moles/mmol * 1e-12 Tg/g * 365 days/year

%%
%March 2017:
%try to calculate "residual" OM in this domain, now instead assuming that
%all o2 coming into the box will be used.

%POCbdiff = 
O2bdiff = O2in.*phib.*vol; %this is the total O2 within phi. from O2div within phi
%O2bdiff = O2div.*phib.*vol; %this is the total O2 within phi. from O2div within phi

residPOC = POCbdiff - O2bdiff*r; %how much "leftover" POC?
%Nloss=residPOC*rN; %mmoles NO3 per day needed to oxidize all that POC
%Nloss = nansum(nansum(nansum(Nloss))); %N lost to C oxidation
%Nloss = Nloss*14*1e-3*1e-12*365 %Tg N/year (from mmoles N *14g/mol * 1e-3moles/mmol * 1e-12 Tg/g * 365 days/year

