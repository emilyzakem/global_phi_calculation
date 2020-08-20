function[diffp]=OCCAdiffp_xyPart1(C,kxall,kyall,dxall,dxvall,dy,dzf)

%2/21/17: now, after rewriting the div file, convert BACK to the p by
%uncommenting the conditionals.

%12/5/16: from OCCAdiffp_xyonly, commented out the < or > 0 conditionals:

%2/20/17: redo do include ddy(dx) at N and S faces.
%and this only does the first part : ddx(Kp*ddx(C)) and ddy(Kp*ddy(C))

%diffdiv = mol C /m3 /s
%calculated at cell centers

nx = size(C,3);
ny = size(C,2);
nz = size(C,1); %22

diffp=nan*ones(nz,ny,nx); %nz is cell centered depths 

for k=1:nx
    for j=1:ny
        for h=1:nz
  
 %get kx,ky: (there was a mistake here in the logic, fixed 2/20/16, but not using this part anyways:)          
          if size(kxall)==1
              kx=kxall;
              kx1=kx;
          else
              kx=kxall(h,j,k);
              if k~=nx; kx1=kxall(h,j,k+1); 
                else kx1=kxall(h,j,1); end %periodic
          end
              
          if size(kyall)==1
              ky=kyall;
              ky1=ky;
          else
              ky=kyall(h,j,k);
              if j~=ny; ky1=kyall(h,j+1,k); end
          end
          
 %%%%%%%%%%%%%%%%%%%%%%%%%         
  %!at south face:
            if j==1;
                Fn=0; %from box to south of j
                Fni=0;
            else
                Fn = ky*(C(h,j,k)-C(h,j-1,k))/dy; %dC/dy at south/bottom face, 
                if Fn<0 %if this is negative, there's more in the j-1 box: INFLUX
                    Fni=Fn*dxvall(j)*dzf(h); %conc*m^3/s, divide by vol below
                else Fni=0; end
            end
                
              
  %!at north face:
            if j==ny
                Fs=0; %no southward incoming at north face, from north
                Fsi=0;
            else
                Fs=ky1*(C(h,j+1,k)-C(h,j,k))/dy;%southward at north face of j
                if Fs>0
                    Fsi=Fs*dxvall(j+1)*dzf(h);
                else Fsi=0; end
            end
                
            
  %at west face:
            if k==1 %%eastward incoming at west face, from box to west of k
                Fe=kx*(C(h,j,k)-C(h,j,360))/dxall(j); %periodic: box 360 is east
            else
                Fe=kx*(C(h,j,k)-C(h,j,k-1))/dxall(j);
            end
                if Fe<0 %negative means 360 box is bigger: INFLUX
                    Fei=Fe*dy*dzf(h);
                else Fei=0; end  
  
  %at east face:
            if k==nx %%westward flux from box to east (box 1) so k1==1
                Fw=kx1*(C(h,j,1)-C(h,j,k))/dxall(j); %k is 360, so k+1 = 1st box
            else
                Fw=kx1*(C(h,j,k+1)-C(h,j,k))/dxall(j); %westward at east face of k
            end            
                if Fw>0 %negative means 360 box is bigger: INFLUX
                    Fwi=Fw*dy*dzf(h);
                else Fwi=0; end 
         
         vol=dxall(j)*dy*dzf(h);
                 
         diffp(h,j,k)=(Fsi-Fni+Fwi-Fei)/vol; %negative fluxes are subtracted here.

          end
        end
    end


end