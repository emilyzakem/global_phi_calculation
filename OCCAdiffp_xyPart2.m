function[diffp]=OCCAdiffp_xyPart2(C,kwx,kwy,dxall,dy,dzf,dzc)
%%

%2/21/17: isolate the GROSS flux
%need to figure out what's incoming and outgoing...

%PArt TWO: ddx(kmx*ddz(C)) and ddy(kmy*ddz(C))

%6/14/20: following xyPart1, add conditionals

%debug:
% kwx=gmkwx;
% kwy=gmkwy;
% dxall=dx;
% dy=dy;
% dzf=dzocca;
% dzc=dzcocca;
% C=O2interp_L;

%i do this with three loops: 1. calc ddz(C); 2. take ddx and ddy at top faces; 3.
%interpolate to horizontal center

%kwx and kwy are defined at cell TOP FACES, i assume, following w and psi,
%which gives me no adv divergence, and hwich is DIFFERENT than documented
%for OCCA (but makes more sense)

%diffdiv = mol C /m3 /s
%calculated at cell centers

nx = size(C,3);
ny = size(C,2);
nz = size(C,1); %22

ddzC=nan*ones(nz,ny,nx); %nz is cell centered depths 
%divFx=nan*ones(nz,ny,nx); %nz is cell centered depths 
%divFy=nan*ones(nz,ny,nx); %nz is cell centered depths 
inFx=nan*ones(nz,ny,nx); %nz is cell centered depths 
inFy=nan*ones(nz,ny,nx); %nz is cell centered depths 

diffp=nan*ones(nz,ny,nx); %nz is cell centered depths 


for k=1:nx
    for j=1:ny
        for h=2:nz
            
 %first, calculate ddz of tracer, which puts it at the same TOP face as kwx and kwy 

        ddzC(h,j,k) = ( C(h,j,k) - C(h-1,j,k) )/dzc(h-1); %dzc at h-1 since we are talking about the TOP face (undefined dzc for h=1) 
        
        end
    end
end

%Jun 2020:
%EJZ: here, the INCOMING fluxes are when ddzC < 0, when h-1 has larger conc (bc ddzC > 0 if conc of box h is higher)

%Redo: label + and n fluxes here (NOT incoming and outgoing, but opposite of that):
ddzCp = ddzC;
ddzCp(ddzC < 0) = 0;

ddzCn = ddzC;
ddzCn(ddzC > 0) = 0;

 %multiply this by kwx(h,j,k) and kwy(h,j,k), defined at same loc as ddzC, and then take x and y gradients:
 %calc ddx and ddy as the CENTERED difference. this gives the flux at the center for that box.
 %then, interpolate from top and bottom to horizontal centers.
 
 %centered difference at TOP FACE of kw*ddzC in x and y: account for differences in dx
 %in dyflux here, THEN divide by vol after:
 %assume that i use centered kwx:
 
 %2/21/17: to get GROSS flux, isolate the POSITIVE contributions 
 %6/24/20: i think the best way to do it is above, since below is the
 %divergence. this is consistent with xyPart1, where outgoing fluxes were
 %zeroed out and then divergence was calculated

 for k=1:nx
    for j=1:ny
        for h=2:nz

%             if k==1
%                 %divFx(h,j,k) =  kwx(h,j,k) * ( ddzC(h,j,k+1) - ddzC(h,j,nx) )*dy*dzc(h-1);
%                 inFx(h,j,k) =  kwx(h,j,k) * ( ddzCn(h,j,k+1) + ddzCn(h,j,nx) )*dy*dzc(h-1);
%             elseif k==nx
%                 %divFx(h,j,k) =  kwx(h,j,k) * ( ddzC(h,j,1) - ddzC(h,j,k-1) )*dy*dzc(h-1);
%                 inFx(h,j,k) =  kwx(h,j,k) * ( ddzCn(h,j,1) - ddzCp(h,j,k-1) )*dy*dzc(h-1);
%             else
%                 %divFx(h,j,k) =  kwx(h,j,k) * ( ddzC(h,j,k+1) - ddzC(h,j,k-1) )*dy*dzc(h-1);
%                 %EJZ Jun 2020: collect incoming based on lat diff
%                 inFx(h,j,k) =  kwx(h,j,k) * ( ddzCn(h,j,k+1) - ddzCp(h,j,k-1) )*dy*dzc(h-1);
%             end

%THIS is a more straightforward way to do it.
%NOTE that kwx and kwy can be + and - !!! SO, then, ignore the negative
%result:
            if kwx(h,j,k) > 0
                inFx(h,j,k) =  kwx(h,j,k) * ( -ddzCn(h,j,k) )*dy*dzc(h-1);
            else
                inFx(h,j,k) =  kwx(h,j,k) * ( -ddzCp(h,j,k) )*dy*dzc(h-1);
            end
            

%             if j==1
%                 %divFy(h,j,k) =  kwy(h,j,k) * ( ddzC(h,j+1,k) - 0 )*dxall(j)*dzc(h-1);
%                 inFy(h,j,k) =  kwy(h,j,k) * ( ddzCn(h,j+1,k) - 0 )*dxall(j)*dzc(h-1);
%             elseif j==ny
%                 %divFy(h,j,k) =  kwy(h,j,k) * ( 0 - ddzC(h,j-1,k) )*dxall(j)*dzc(h-1);
%                 inFy(h,j,k) =  kwy(h,j,k) * ( 0 - ddzCp(h,j-1,k) )*dxall(j)*dzc(h-1);
%             else
%                 %use centered dx since i am calculating it at the center... (??)
%                 %divFy(h,j,k) =  kwy(h,j,k) * ( ddzC(h,j+1,k) - ddzC(h,j-1,k) )*dxall(j)*dzc(h-1);
%                 %EJZ Jun 2020: collect incoming based on lat diff
%                 inFy(h,j,k) =  kwy(h,j,k) * ( ddzCn(h,j+1,k) - ddzCp(h,j-1,k) )*dxall(j)*dzc(h-1);
%             end

            if kwy(h,j,k) > 0
                inFy(h,j,k) =  kwy(h,j,k) * ( -ddzCn(h,j,k) )*dxall(j)*dzc(h-1);
            else
                inFy(h,j,k) =  kwy(h,j,k) * ( -ddzCp(h,j,k) )*dxall(j)*dzc(h-1);
            end
            
        end
    end
 end
 
  %now interpolate from top face to horizontal center (h+1 + h)/2, which works since center is
 %center

 for k=1:nx
    for j=1:ny
        for h=1:nz-1
            
       %divFxc = (divFx(h+1,j,k) + divFx(h,j,k))/2;
       %divFyc = (divFy(h+1,j,k) + divFy(h,j,k))/2;
       inFxc = (inFx(h+1,j,k) + inFx(h,j,k))/2;
       inFyc = (inFy(h+1,j,k) + inFy(h,j,k))/2;

    
 %then, the total div from these two terms is the divergence in x and the
 %div in y. here, divide by volume to complete the ddx and ddy portions:
 
        vol=dxall(j)*dy*dzf(h);
 
        %diffdiv2(h,j,k) = (divFxc + divFyc)/vol;
        diffp(h,j,k) = (inFxc + inFyc)/vol;
 
        
        end
    end
 end
%%

end