function[diffp_z2]=ddz_kz_expl2_p(C,kz2,dxall,dy,dzf)
%EJZ: on 6/25/20, try resolving the KPP and GMkwz explicitly, rather than
%relying on the implicit solver, to check for dependency on that solver.

nx = size(C,3);
ny = size(C,2);
nz = size(C,1); %22

kz2_explicit=nan*ones(nz,ny,nx);
diffp_z2=nan*ones(nz,ny,nx);

%calc flux at top face
kz2_explicit(2:end,:,:) = kz2(2:end,:,:).*( C(2:end,:,:) - C(1:end-1,:,:) )./dzf(2:end,:,:);

df = kz2_explicit;
dfp=df;
dfn=df;

dfp(dfp<0)=0; %positive fluxes IN at bottom face (h+1)
dfn(dfn>0)=0; %negative fluxes IN at top face (h)

diffp_z2(1:end-1,:,:) = (dfp(2:end,:,:) - dfn(1:end-1,:,:))./dzf(1:end-1,:,:);

end


