%6/12/20: modified (annotated) from O2inNew_Load.m by EJZ

%2/21/17:
%from the div script O2inNew_Load, convert this to the INCOMING gross flux
%of oxygen. 

%ECCO velocity data

%OCCA v 2, annual resolution: 50 depth levels: 50x160x360
%http://mit.ecco-group.org/opendap/ecco_for_las/
%annual
ncload('OCCAvelocity/DDuvel.0406annclim.nc.nc') %depth_c x latitude_t x long_u  (size: 23x160x360)
uecco=u; %m/s
%long u is at WESTERN edge of box (0:359) while long t is 0.5 to 359.5

ncload('OCCAvelocity/DDvvel.0406annclim.nc.nc') %depth_c x latitude_v x long_t   (size: 23x160x360)
vecco=v; %m/s
%v is at SOUTHERN edge of box (-80 to 79), whereas lat_t is -79.5 to 79.5

ncload('OCCAvelocity/DDwvel.0406annclim.nc.nc')
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

%deconstruct: (the sum of this is not exactly advp) XX delete this if this
%is OK.
advp_mean=OCCAadvp(O2interp_L,uecco,vecco,wecco,dx,dxv,dy,dzocca);
advp_bolus=OCCAadvp(O2interp_L,ubg,vbg,wbg,dx,dxv,dy,dzocca);

%ones:
advp_ONES=OCCAadvp(ONESi,uadd,vadd,wadd,dx,dxv,dy,dzocca);

%%
%Horizontal Mixing:

disp('Diffusion in x and y...')

%Part 1 (simple 1e3)
diffp_xy1=OCCAdiffp_xyPart1(O2interp_L,1e3,1e3,dx,dxv,dy,dzocca); %_xyPart1

%Part 2 (GM)
diffp_xy2=OCCAdiffp_xyPart2(O2interp_L,gmkwx,gmkwy,dx,dy,dzocca,dzcocca);

%%

%Vertical Mixing:

disp('Vertical mixing...')

%Part 1 (explicit from GM)

%this calculates: ddz( [kwx]*(dC/dx) + [kwy]*(dC/dy) )
%diffp_z1=ddz_kz_expl1(O2interp_L,gmkwx,gmkwy,dx,dy,dzocca);
diffp_z1=ddz_kz_expl1_p(O2interp_L,gmkwx,gmkwy,dx,dy,dzocca); %6/24/20


%Part 2 (implicit from KPP and gmkwz)

disp('Implicit vertical mixing...')

kz_2=kpp+gmkwz; %kpp part here?? check with Jon

%turn nans to zeros:
kz_2(isnan(kz_2))=0;

O2interp0=O2interp_L;
O2interp0(isnan(O2interp0))=0;

%first, new dimensions for the mex-fortran files:
% dims in order (for vars, too) are: integer*4 iMax,jMax,Nr,Nt
dims=[size(O2interp_L,3) size(O2interp_L,2) size(O2interp_L,1) 1];

O2interp0f=permute(O2interp0,[3 2 1]);
O2interp0f2=repmat(O2interp0f,[1 1 1 2]);
O2interp0f=O2interp0f2(:,:,:,1);

kz_2f=permute(kz_2,[3 2 1]);
kz_2f2=repmat(kz_2f,[1 1 1 2]);
kz_2f=kz_2f2(:,:,:,1);

%make grids of ones for the mex file
gridones=ones(dims(1:3));
surfones=ones(dims(1:2));

%get dzc as nz long:

%dzc is the z between cell centers, starting with the diff between box 1
%and box 2 (k and k+1)
dzc1=[0; dzcocca]; %so i have k centered dzs - this is the right end to put nan on bc the code calls for dzc(k)
%dzc(k) for C(k)-C(k-1) and dzc(k+1) for C(k+1)-C(k) face

%Calculate implicit vertical mixing using  MEX
%output is the FLUX in z: wtidiff = df = mol/m2/s
%timestep=43200; % 12hr timestep. from john.
%try lower timestep: (one hour, why not:) %
timestep=3600;
%iterates 12 times, possibly should come up with a criteria to stop
%this is the implicitly calculated vertical diffusive flux: mol O2 /m2 /s at each face.
%is it at the top or bottom face? code seems like it's the flux between k
%and k-1, so, that would imply that this is the flux at the TOP face:
%2/20/17: i agree this should be the flux at the top face. this agrees with
%the way psi and w and gmkwz are defined (which is different than OCCA
%documentation)

[df,xx]=mit_difimpl_flux(timestep,surfones,dzocca,dzc1,gridones,...
    kz_2f2,O2interp0f2); %mol/s which is the variable 'df,' the diffusive fluxes
%df is now 360x160x50x2

%2/20/17: df is the diffusive FLUX. it is zero at the top face.

%this means df(1) should be zero (no flux from top into first box)
%and df(end) is the flux into the last box from the second to last box.
%is positive or negative in or out? i think positive is UP.
%if so, flux at k > 0 goes into box k-1. flux at k < 0 goes into box k

%i just want the last time step (the second of the two):
%i ASSUME this. but not really sure what these two timesteps mean, if
%anything.
df = df(:,:,:,2); %diffusive flux at TOP face.
df(df==0)=nan;
dfp=df;
dfn=df;

dfp(dfp<0)=0; %positive fluxes
dfn(dfn>0)=0;

dfpu=zeros(size(dfp));
dfpu(:,:,1:end-1)=dfp(:,:,2:end);

%make dz grid:
dzfgrid = repmat(dzocca,1,dims(1),dims(2));
dzfgrid = permute(dzfgrid,[2 3 1]);

%This gives only the positive fluxes in:
diffp_z2 = (dfpu - dfn)./dzfgrid; 

%resize
diffp_z2 = permute(diffp_z2,[3 2 1]);

%%

%ALTERNATIVELY, try explicit instead of implicit for z2:
%kz_22=kpp+gmkwz;
%diffp_z2=ddz_kz_expl2_p(O2interp_L,kz_22,dx,dy,dzocca); %6/24/20



%%

%cumulative impact of diff terms on N loss:
%with just xy1: 57.5 (old)
%with xy1 + xy2: (with incoming ddzC isolated 54.2), with new way: 54.1 !!
%consistent
%with xy1 + z2: 46
%with xy1 + xy2 + z2: 45.1
%with xy1 + z1: 53
%all four: 43.3!!! phiU: 17, phiD: 103

%diffp = diffp_xy1 + diffp_z2;%+ diffp_xy2 + diffp_z2;% + diffp_z;

%diffp = diffp_xy1 + diffp_xy2 + diffp_z2;

diffp = diffp_xy1 + diffp_xy2 + diffp_z1 + diffp_z2;

%%
        
O2advp=-advp;  %flux units: mmol/m3/s = m/s * mmol/m3 * 1/m
O2diffp=diffp;

%change units: from mmol/m3/s:
O2advp=O2advp*3600*24; %mmol/m3/day = uM O2/day
O2diffp=O2diffp*3600*24; %mmol/m3/day = uM O2/day
O2in_L=O2advp+O2diffp; %mmol/m3/day = uM O2/day

ONESadvp=-advp_ONES*3600*24; %1/day %12/7 %obv no diff for ones bc homogenous (checked: it is 0)
ONESin_L=ONESadvp; 