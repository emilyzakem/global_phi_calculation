%12/8/16: ejz: modeled more JUST off of gmtransport, after comparing with
%Jon Lauderdale's. also, big addition: i divide by grid faces below to
%actually give velocities in m/s as the output. 

%modeled off of mit_gmtransport AND mit_gmvel:
function [ubg,vbg,wbg]=ejz_gmvel(isAdvForm,gmpsix,gmpsiy,dxall,dxvall,dy,dzf)
% Caluculate MITgcm Bolus advective transport for the GM-eddy mixing
% parameterisation. Kwx,Kwy and Kwz are Kgm*Isopycnal_Slope (skew flux) or
% use GM_PSIX and GM_PSIY for the advective form.

%12/8/16: gmpsix and y are in m2/s
%2/18/16: and are on the LOWER LEVEL of the grid: zw (Depth_l) according to
%netcdf outputs        
%2/19/16: psi is NOT on the lower level of the gridbox even though it SAYS
%it is in the documentation. same for w. divergence goes to zero if i
%consider them on the UPPER level of the grid box!!
%also, need to consider dxv -- the difference in dx at the faces vs the
%middle of the box

%for debugging:
%isAdvForm=1;
%dxall=dx; %calc at latt, not latv!
%dzf=dzocca; %depths of each grid box (50 long)

nx = size(gmpsix,3);
ny = size(gmpsix,2);                             
nz = size(gmpsix,1); %22

ubg=nan*ones(nz,ny,nx);
vbg=nan*ones(nz,ny,nx);

%gmadvtest=mit_getparm(data_gmredi,'GM_AdvForm');

%if strcmpi(gmadvtest,'true') % Using advective form
if isAdvForm==1;
% Setting of KGM. Kwx and Kwy are the slope*thickness diffusivity so why are we multiplying by KGM again?
    % This is necessary because the gmredi tensor when the skew flux is used and Kgm=Krho is:
    %    x,    y,   z
    % u: 1,    0,   0
    % v: 0,    1,   0
    % w: 2Sx,  2Sy, |S|^2
    % So to get the slopes to calculate the bolus transport we need to divide by 2.
    % However, the Advective form is used then kwx=Sx and kwy=Sy, hence derived bolus
    % velocities will be half as great as expected.
    kgm=2;
else
    kgm=1;
end

%for k=1:nx
 %   for j=1:ny
       % for h=1:nz %POC flux # boxes in vertical (22)
  %          
            %%From the transport:
            %ubg(i,j,k,t) = (kgm/2).*grid.dyg(i,j)*(GM_PsiX(i,j,kp1,t)*maskp1-GM_PsiX(i,j,k,t))*grid.hfacw(i,j,k);
            %vbg(i,j,k,t) = (kgm/2).*grid.dxg(i,j)*(GM_PsiY(i,j,kp1,t)*maskp1-GM_PsiY(i,j,k,t))*grid.hfacs(i,j,k);
                    
            %%if h~=nz; %ubg IS one larger than this! so code works:
            %ubg(h,j,k) = (kgm/2).*dy*(gmpsix(h+1,j,k)-gmpsix(h,j,k));
            %vbg(h,j,k) = (kgm/2).*dxall(j)*(gmpsix(h+1,j,k)-gmpsix(h,j,k));
            %%end
                            
%             %wbg(i,j,k,t) =(...
%             %                 grid.dyg(ip1,j)*GM_PsiX(ip1,j,k,t)...
%             %                -grid.dyg( i ,j)*GM_PsiX( i ,j,k,t)...
%             %                +grid.dxg(i,j+1)*GM_PsiY(i,j+1,k,t)...
%             %                -grid.dxg(i, j )*GM_PsiY(i, j ,k,t)...
%             %                ).*grid.hfacc(i,j,k);
%             wbg(h,j,k) =  dy*gmpsix(h,j,kp1) ...
%                          -dy*gmpsix(h,j,k) ...
%                          +dxall(j)*gmpsiy(h,j+1,k) ...
%                          +dxall(j)*gmpsiy(h,j,k);
            %%
%my velocities: adapted from mit_gmvel.m above       

%2/19/17: realized that w is actually at the top. and so psi also. so
%documentation may be wrong. 

dxallgr=repmat(dxall,1,nx,1);
dxallgr=permute(dxallgr,[3 1 2]);

    for h=1:nz-1 %psi at top, like w:
        
            ubg(h,:,:) = (kgm/2)*(gmpsix(h+1,:,:)-gmpsix(h,:,:))/dzf(h);
            vbg(h,:,:) = (kgm/2)*(gmpsiy(h+1,:,:)-gmpsiy(h,:,:))/dzf(h);
    end                
            ubg(nz,:,:) = (kgm/2)*(0-gmpsix(h,:,:))/dzf(h);
            vbg(nz,:,:) = (kgm/2)*(0-gmpsiy(h,:,:))/dzf(h);

%vertical velocity:

%part 1:
vbcy=nan*ones(nz,ny,nx);

    for k=1:nx
      for j=1:ny
        for h=1:nz %back to psi at top
            
           if j~=ny %northern edge
                vbcy(h,j,k) = gmpsiy(h,j+1,k)*dxvall(j+1)-gmpsiy(h,j,k)*dxvall(j);
           %else
           %     vbcy(h,j,k) = 0-gmpsiy(h,j,k)*dxall(j); %check this- in code looks like it'd not calculated
           end

        end
      end
    end
           
%part 2:    
ubcx=nan*ones(nz,ny,nx);

    for k=1:nx
      for j=1:ny
        for h=1:nz
            
                if k~=nx
                    kp1=k+1;
                else
                    kp1=1; %periodic in e-w
                end
                ubcx(h,j,k) = (gmpsix(h,j,kp1)-gmpsix(h,j,k))*dy;
       
        end
      end
    end           
     
%%
wbg=vbcy+ubcx;
            
%now divide by faces to get velocities: from m3/s to m/s:

dxallgr2=repmat(dxallgr,nz,1,1);
%dxallgr2=repmat(dzf,1,ny,nx);

wbg = wbg./dy./dxallgr2;

end
  