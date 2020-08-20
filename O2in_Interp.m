%2/21/17, modified from div file

%interpolate O2in and POCdiv onto the same grid:
%POC flux is ALREADY interpolated to the same horizontal grid!

%need to run the two above scripts first to get the right grid variables
%for this script
latt; %same as latosh: -79.5:1:79.5 (centered)
lont; %same as lono: 0.5:1:359.5 (centered)

%POC grid:
zf=unique(pdata(:,3)); %depth at TOP FACES of boxes 
dzf;
zc;

%%
%now get to oxygen incoming:
%INTERPOLATE:

%first, get velocity DEPTH res to same size as POC flux depth res: 
    %i should put velocities at the CENTERS.
    O2in=nan*ones(length(zc),length(latt),length(lont)); %at zcentered
    ONESin=nan*ones(length(zc),length(latt),length(lont)); %at zcentered
    O2interp=nan*ones(length(zc),length(latt),length(lont)); %at zcentered

    disp('Interpolating...')
    for k=1:360
        for j=1:160
            
            kj=squeeze(O2in_L(:,j,k));
            F=griddedInterpolant(zt,kj);
            O2in(:,j,k)=F(zc);
            
            kj=squeeze(ONESin_L(:,j,k));
            F=griddedInterpolant(zt,kj);
            ONESin(:,j,k)=F(zc);
            
            kj=squeeze(O2interp_L(:,j,k));
            F=griddedInterpolant(zt,kj);
            O2interp(:,j,k)=F(zc);
            
        end  
    end
    disp('Done interpolating.')
    
    %%
    
    figure;
    set(gcf,'color','w')
    
    subplot(1,2,1)
    pcolor(lont,latt,squeeze(O2in(4,:,:)))
    shading flat
    caxis([-1 1])
    colorbar
    makerbcolormap;
    title('O2in at 270 m')
    set(gca,'fontsize',16)

    subplot(1,2,2)
    pcolor(lont,latt,squeeze(POCdiv(4,:,:)))
    shading flat
    caxis([-.1 .1])
    colorbar
    makerbcolormap;
    title('POCdiv at 270 m')
        set(gca,'fontsize',16)
