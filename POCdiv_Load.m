

%%
%POC FLUX data
pdata=importdata('POC_FLUXES_Schlitzer-2002_Exp_C_DBI1niA.txt');
pdata=pdata.data; %lon (degE), lat (degN), depth (m), POC flux (mol/m2/year)

%find unique values for discrete points:
lonp=unique(pdata(:,1)); %non-uniform
latp=unique(pdata(:,2)); %non-uniform -MOST are at less often points than this. (ex: -12, not -11)

%DEPTH:
zf=unique(pdata(:,3)); %depth at TOP FACES of boxes (where POC flux is defined). uniform here!
%dz at faces (height of boxes):
dzf=zf(2:end)-zf(1:end-1); %only 22 long here... ignore deepest box in long run bc of this.. (23rd face of POC flux)
zbot=zf(2:end); %bottom depth (since i'm essentially ignoring the box below zf(end))

%Centered depth:
zc=zf(1:end-1)+dzf/2; %this is also only 22 long
%dz at centers:
dzc=(dzf(1:end-1)+dzf(2:end))/2; %only 21 long here, bc we're not dealing with top or bottom fluxes.
    
%NOW, I have to interpolate here bc i can't figure out how else to do it...

poc3d=nan*ones(length(zf),length(latp),length(lonp));
poc3di=nan*ones(length(zf),length(latt),length(lont));

 %to make land mask, interpolate O2 by depth TEMPORARILY to velocity/POC centered depths:
    O2interpmask=nan*ones(length(zc),size(uecco,2),size(uecco,3));

    for k=1:360
        for j=1:160
            o2kj=squeeze(O2sh(:,j,k));
            F=griddedInterpolant(zo,o2kj);
            O2interpmask(:,j,k)=F(zc);
        end
    end

    %%
for h=1:length(zf)-1 %don't put values on very last level

list=pdata(pdata(:,3)==zf(h),:);

%%
pocg=nan*ones(length(latp),length(lonp));

for j=1:length(lonp);
    for k=1:length(latp);        
        lonj=lonp(j);
        latk=latp(k);

        pjk=list(list(:,1)==lonj&list(:,2)==latk,4);
        if isempty(pjk)
            pjk=nan;
        end        
        pocg(k,j)=pjk;
    end
end
%%
%interpolate:

F=scatteredInterpolant(list(:,1),list(:,2),list(:,4));

[xg,yg]=meshgrid(lono,latosh);
%nan out the land values:
if h<23
omask=squeeze(O2interpmask(h,:,:));
else
    omask=squeeze(O2interpmask(22,:,:));
end

omask(~isnan(omask))=1;
yg=yg.*omask;
xg=xg.*omask;

p2=F(xg,yg);

poc3d(h,:,:)=pocg;
poc3di(h,:,:)=p2; %mol C/m2/year

end

%%
%plot to match henson 2011:
figure;
plotd=5;%3 for 133m, 14 for 2060m
pcolor(lont,latt,-squeeze(poc3di(plotd,:,:))*12); 
title(['Schlitzer POC flux: ' num2str(zf(plotd)) 'm (gC/m2/yr)'],'fontsize',16)
%pcolor(lono,latosh,-squeeze(poc3di(3,:,:))*12*1e3/365);
%title('Schlitzer POC flux (mgC/m2/d)')
shading flat
colorbar
colormap jet
          %colormap jet, cmap=colormap; colormap(cmap([15:3:59 62],:));

%caxis([0 8])
set(gcf,'color','white')
set(gca,'fontsize',14)



%%
%end???
%%
%now calculate the incoming POC flux in mol C/m3/year by dividing by
%length:
dzfgrid = repmat(dzf,1,size(poc3di,2),size(poc3di,3));
%dzfgrid = permute(dzfgrid,[2 3 1]);

%12/5/16:
%divide by depth and take out bottom layer (all nans anyways):
%just in: POCin=poc3di(1:end-1,:,:)./dzfgrid; %mol C/m3/year
%POCin=POCin*1e3/365; %from molC/m3/year to mmolC/m3/day = uM/day

%to get divergence of POC, take difference:
%from before: POCbdiff=POCin(2:end,:,:)*stuff-POCin(1:end-1,:,:)*stuff; %left in each box
POCdiv = (poc3di(2:end,:,:) - poc3di(1:end-1,:,:))./dzfgrid; %POSITIVE: no POC gained
POCdiv=POCdiv*1e3/365; %from molC/m3/year to mmolC/m3/day = uM/day

