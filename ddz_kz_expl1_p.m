function[diffp_z1]=ddz_kz_expl1_p(C,kwx,kwy,dxall,dy,dzf)
%modified from:
%function[kz_explicit]=kz_expl1(C,kwx,kwy,dxall,dy)
%on 2/19/17
%calculate [kwx]*(dC/dx) and [kwy]*(dC/dy)
%and then take the vertical gradient of it: ddz( [kwx]*(dC/dx) + [kwy]*(dC/dy)  )

%2/19/17 NOTE: I need to check that the SIGN is right for this: that ddz is
%positive or negative in the right direction....

%Jun 2020: _p means positive fluxes are taken now

%this is the factor: [kwx]*(dC/dx) = Kp*Sx*(dC/dx)

nx = size(C,3);
ny = size(C,2);
nz = size(C,1); %22

%for centered differences of tracer:
Cdiffx=nan*ones(nz,ny,nx); %2/19/17: ideally this should be calculated at the TOP face
Cdiffy=nan*ones(nz,ny,nx); %2/19/17: ideally this should be calculated at the TOP face

kz_explicit=nan*ones(nz,ny,nx); %2/19/17: ideally this should be calculated at the TOP face
%assuming that kmx,kmy,kmz are at TOP depth, not Depth_l, and that
%documentation is wrong: same as w and GMpsix and GMpsiy

diffp_z1=nan*ones(nz,ny,nx);

for k=1:nx
    for j=1:ny
        for h=2:nz %POC flux # boxes in vertical (22)
            
            if j~=ny %none coming out of top of northern edge
                if j~=1 %or southern
                    if k~=nx
                        kp1=k+1;
                    else
                        kp1=1; %periodic in e-w
                    end
                    
                    if k~=1
                        km1=k-1;
                    else
                        km1=nx; %periodic in e-w
                    end
                    %C is at cell centers, but kwx and kwy are at cell bottom face in center of that.
                    %so take centered difference of tracer to get the resulting dC/dx at the
                    %center, too.
                    
                    Cdiffx(h,j,k) = ( C(h,j,kp1) - C(h,j,km1) )/(2*dxall(j));
                    Cdiffy(h,j,k) = ( C(h,j+1,k) - C(h,j-1,k) )/(2*dy); %this doesn't incorporate the different
                    %areas of the N and S faces, but, this is just the tracer difference, so accurate
                    
                    %since i ultimately want to get a vertical gradient of this that ends up centered (at h+1/2),
                    %i want to calculate this Cdiffx and Cdiffy at the top and bottom faces:
                    %so get the average at the top edge (h and h-1) for kw at h, and avg at bottom edge (h and h+1) for kw at h+1
                    
                    %calc kz at TOP face
                    kz_explicit(h,j,k) = ...
                        kwx(h,j,k)*( Cdiffx(h,j,k)+Cdiffx(h-1,j,k) )/2 ...
                        +kwy(h,j,k)*( Cdiffy(h,j,k)+Cdiffy(h-1,j,k) )/2;
                    
                    %OLD way, not centering dC/dx and dC/dy
                    %kwx(h,j,k)*( C(h,j,kp1)-C(h,j,k) )/dxall(j) ...
                    %+kwy(h,j,k)*( C(h,j+1,k)-C(h,j,k) )/dy;
                    
                end %j not 1
            end %j not ny
            
        end
    end
end

%now take vertical gradient of this. kz_explicit is defined now at the top
%face of each cell. so gradient needs to be a f of h+1 and h, with dz being
%the height of the box itself (not the height between cell centers) which
%is dzocca (here as dzf)

% for k=2:nx
%     for j=1:ny
%         for h=1:nz-1 %POC flux # boxes in vertical (22)
%
%             %is the sign correct here? CHECK against other diffusion routine...
%             %2/19/17 checked: if Cdiffx is POS, means C outside of box is larger, so flux in.
%             %then, if that is POS at h+1 (and smaller at h), then the below has a flux
%             %going INTO the box. if this is POS at h (and smaller at h+1), this is a
%             %sink from the box: the dz gradient is <0. i think this SIGN is right.
%
%             %OLD DIVERGENCE:
%             %diffp_z1(h,j,k) = ( kz_explicit(h+1,j,k) - kz_explicit(h,j,k) )/dzf(h);
%
%             %EJZ Jun 2020: which one is incoming? assume this is the same as the
%             %implicit, and that a positive df (kz_explicit > 0) means flux is going UP,
%             %or out of box. following implicit step:
%
%
%         end
%     end
% end

%EJZ Jun 2020: which one is incoming? assume this is the same as the
%implicit, and that a positive df (kz_explicit > 0) means flux is going UP,
%or out of box. following implicit step:



df = kz_explicit;
dfp=df;
dfn=df;

dfp(dfp<0)=0; %positive fluxes IN at bottom face (h+1)
dfn(dfn>0)=0; %negative fluxes IN at top face (h)

diffp_z1(1:end-1,:,:) = (dfp(2:end,:,:) - dfn(1:end-1,:,:))./dzf(1:end-1,:,:);



end


