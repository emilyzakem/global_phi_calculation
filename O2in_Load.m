%6/12/20: modified (annotated) from O2inNew_Load.m by EJZ

%2/21/17:
%from the div script O2inNew_Load, convert this to the INCOMING gross flux
%of oxygen. 

%ECCO velocity data

%OCCA v 2, annual resolution: 50 depth levels: 50x160x360
%http://mit.ecco-group.org/opendap/ecco_for_las/
%annual
ncload('/Users/ezakem/Documents/DeepRespiration/data/O2andECCOandPOCflux/OCCAvelocity/DDuvel.0406annclim.nc.nc') %depth_c x latitude_t x long_u  (size: 23x160x360)
uecco=u; %m/s
%long u is at WESTERN edge of box (0:359) while long t is 0.5 to 359.5

ncload('/Users/ezakem/Documents/DeepRespiration/data/O2andECCOandPOCflux/OCCAvelocity/DDvvel.0406annclim.nc.nc') %depth_c x latitude_v x long_t   (size: 23x160x360)
vecco=v; %m/s
%v is at SOUTHERN edge of box (-80 to 79), whereas lat_t is -79.5 to 79.5

ncload('/Users/ezakem/Documents/DeepRespiration/data/O2andECCOandPOCflux/OCCAvelocity/DDwvel.0406annclim.nc.nc')
wecco=w; %m/s  depth_l x latitude_t x long_t  (size: 50x160x360)
%w is at BOTTOM depth of grid cell (10 to 6134.5) whereas depth_c is 5 to
%5906

%nan out missing values:
uecco(uecco<-1e20)=nan;
vecco(vecco<-1e20)=nan;
wecco(wecco<-1e20)=nan;

latv=Latitude_v; % -80:1:79 ("south" face)
latt=Latitude_t; % -79.5:1:79.5 (centered)
lonu=Longitude_u; % 0:1:359 ("west" face)
lont=Longitude_t; % 0.5:1:359.5 (centered)
%zw=Depth_w; % 0:irregular:5700 51 points, all edges)
zw=Depth_l; %BOTTOM depths: 50 long. ends up with new ecco (occa): w does not have one extra edge!!!!!
zt=Depth_c;  % 5:irregular:5450 (50 points, "centered")

%lengths of depths:
dzocca=Depth_w(2:end)-Depth_w(1:end-1); %between faces (50 long)
dzcocca=(dzocca(1:end-1)+dzocca(2:end))/2; %between centers (49 long)


%velocity is in m/s
%dx is a f of Y, not x, and so, dx at faces = dx at centers:
dx0 = 6.371e6*pi/180;
dx = dx0*cos(latt*pi/180);
%dx=111.32*1e3*cos(latt*pi/180); %m width of 1 degree = 

%dx at south side of grid cell, for gmpsi calc;
dxv=dx0*cos(latv*pi/180); %m width of 1 degree = 

%dy=111.32*1e3; %latitudinal distance is relatively constant from N to S
dy=dx0;
%%
%MIXING:
%OCCA paper: Krho is 1000
%Vertical diffusion coefficient for salt & tracer:
%m2/s latt,lont, zt (depth_c)
ncload('OCCAvelocity/DKPPdifs.0406annclim.nc.nc')
kpp=kppdifs; %m^2/s
kpp(kpp<-1e20)=nan;
%kpp=abs(kpp); %no neg values- this was ON to run the 3.14.16 diff values.

%%
%GM:
%from Jon Lauderdale (3/10/16): mit_gmtransport.m
%psix: m2/s: lonu, latt, zw (depth_l)
ncload('OCCAvelocity/DGMPsiX.0406annclim.nc.nc')
%psiy: m2/s: lont, latv, zw (depth_l)
ncload('OCCAvelocity/DGMPsiY.0406annclim.nc.nc')
gmpsix(gmpsix<-1e20)=nan;
gmpsiy(gmpsiy<-1e20)=nan;
%%
%kwx,y,z: m2/s with lont, latt, zw (depth_l at bottom of box):
%kwz is positive definite. kwx and kwy are NOT:
ncload('OCCAvelocity/DGMkwx.0406annclim.nc.nc')
ncload('OCCAvelocity/DGMkwy.0406annclim.nc.nc')
ncload('OCCAvelocity/DGMkwz.0406annclim.nc.nc')
gmkwx(gmkwx<-1e20)=nan;
gmkwy(gmkwy<-1e20)=nan;
gmkwz(gmkwz<-1e20)=nan;

%%
%OXYGEN:
%new oxygen data from WOA 2013 (updated with higher vertical res):
%WITHOUT Bianchi correction
ncload('woa13_all_o00_01.nc')

lato=lat;
%lono=lon;
zo=depth; %102 long

%rearrange to match velocity:
lono=[lon(181:end); lon(1:180)+180*2];
O2=cat(3,o_an(:,:,181:end),o_an(:,:,1:180));
%O2=permute(O2,[3 2 1]);
O2(O2>9e36)=nan;
O2=O2/22391.6*1000.0*1e3; %from ml/l to uM
%Bianchi correction:
O2=1.009*O2-2.523;
%disp('No correction for O2 ')

%%
estO2error

%%

%INTERPOLATE OXYGEN: get O2 to same depth res as ECCO:
    %and N and P also from WOA 2013
    %first, cut down to same 160 lat band as ecco:
    O2sh=O2(:,11:170,:);   
    latosh=lato(11:170); %this is the lato that corresponds to O2sh and O2in. latosh=latt

%now, get O2 to the same depths as OCCA: we want O2 on centered depths: zt
    %(Depth_c)
        O2interp_L=nan*ones(length(zt),size(uecco,2),size(uecco,3));

    for k=1:360
        for j=1:160
            
            o2kj=squeeze(O2sh(:,j,k));
            F=griddedInterpolant(zo,o2kj);
            O2interp_L(:,j,k)=F(zt);
           
        end
    end
    
    ONESi=wecco;
    ONESi(~isnan(wecco))=1;
    %ONESi(~isnan(PO4i))=1000; %test, 4/11/16   
    
%%
%Now, with O2 interpolated to velocities, we can calculate the oxygen
%divergence.

%get total velocities:

%get bolus transport from psix and psiy:
%modeled off of mit_gmtransport.m:
%function [ubg,vbg,wbg]=ejz_gmvel(isAdvForm,gmpsix,gmpsiy,dxall,dy,dzf)
%psi is on bottom face:
%dx is used for v calc, which is on south face, so use dxv here:
%as of Feb 19, 2017, this is working correctly with no resulting divergence
[ubg,vbg,wbg]=ejz_gmvel(1,gmpsix,gmpsiy,dx,dxv,dy,dzocca); 


%now add bolus to geostrophic velocities:
uadd=uecco+ubg;
vadd=vecco+vbg;
wadd=wecco+wbg;

%now advect:
advp=OCCAadvp(O2interp_L,uadd,vadd,wadd,dx,dxv,dy,dzocca);

%deconstruct:
advp_mean=OCCAadvp(O2interp_L,uecco,vecco,wecco,dx,dxv,dy,dzocca);
advp_bolus=OCCAadvp(O2interp_L,ubg,vbg,wbg,dx,dxv,dy,dzocca);

%ones:
advp_ONES=OCCAadvp(ONESi,uadd,vadd,wadd,dx,dxv,dy,dzocca);

%%
%Horizontal Mixing:

disp('Diffusion in x and y...')

diffp_xy1=OCCAdiffp_xyPart1(O2interp_L,1e3,1e3,dx,dxv,dy,dzocca); %_xyPart1

diffp_xy1_ONES=OCCAdiffp_xyPart1(ONESi,1e3,1e3,dx,dxv,dy,dzocca);


diffp_xy2=OCCAdiffp_xyPart2(O2interp_L,gmkwx,gmkwy,dx,dy,dzocca,dzcocca);

diffp_xy2_ONES=OCCAdiffp_xyPart2(ONESi,gmkwx,gmkwy,dx,dy,dzocca,dzcocca);

%nanbox = nan*ones(size(diffp_xy1));
%diffp_xy2n = nanbox;
%diffp_xy2n(2:end,:,:) = diffp_xy2;
%diffp_xy2_ONESn = nanbox;
%diffp_xy2_ONESn(2:end,:,:) = diffp_xy2_ONES;


%%

%Vertical Mixing:
%compile the kz mixing term:

disp('Vertical mixing...')

%this calculates: ddz( [kwx]*(dC/dx) + [kwy]*(dC/dy) )
diffp_z=ddz_kz_expl1(O2interp_L,gmkwx,gmkwy,dx,dy,dzocca);

diffp_z_ONES=ddz_kz_expl1(ONESi,gmkwx,gmkwy,dx,dy,dzocca);

%%

diffp = diffp_xy1 + diffp_xy2 + diffp_z;

diffp_ONES = diffp_xy1_ONES + diffp_xy2_ONES + diffp_z_ONES; 

%%
        
O2advp=-advp;  %flux units: mmol/m3/s = m/s * mmol/m3 * 1/m
O2diffp=diffp;

%change units: from mmol/m3/s:
O2advp=O2advp*3600*24; %mmol/m3/day = uM O2/day
O2diffp=O2diffp*3600*24; %mmol/m3/day = uM O2/day
O2in_L=O2advp+O2diffp; %mmol/m3/day = uM O2/day

ONESadvp=-advp_ONES*3600*24; %1/day %12/7
ONESdiffp=diffp_ONES*3600*24;
ONESin_L=ONESadvp+ONESdiffp;
