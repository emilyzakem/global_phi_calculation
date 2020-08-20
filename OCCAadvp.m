function[advp]=OCCAadvp(C,u,v,w,dxall,dxvall,dy,dzf)

%12/5/16: found a typo:

%calculates the INCOMING fluxes due to advection of C by u,v,w
%for grid size nk x nj x nh and cell size dxallall(j) x dy x dzall(h)

%2/19/17: fix with the new dxvall (the diff in the dx across the box
%matters) and also use centered tracer concentration.

nx = size(C,3);
ny = size(C,2);
nz = size(C,1); %22

advp=nan*ones(nz,ny,nx);

for k=1:nx
    for j=1:ny
        for h=1:nz %POC flux # boxes in vertical (22)
    
        %advective:
        %adopted from "myupwind" scheme from JunG6_upw5
        
        % j loop: meridional velocity v- starts at south (-79.5 centered latt(1))
                %along latitude- south and north faces- v and j
                %only need south face positive and north face negative velocities:
        
       %Fluxes: v*O2 = m/s * mmol/m3 = mmol/m2/s
            %!at south face:
            if j==1;
                Fn=0; %northward incoming at south face, from box to south of j
            else
                vp=(v(h,j,k)+abs(v(h,j,k)))/2; %for lat= -80 to 79
                %vn=(v(h,j,k)-abs(v(h,j,k)))/2;
                %UPWIND:
                %Fn=vp*C(h,j-1,k); %northward incoming at south face, from box to south of j
                %CENTERED:
                Cf=(C(h,j-1,k)+C(h,j,k))/2;
                Fn=vp*Cf*dxvall(j)*dzf(h);
            end

            %!at north face:
            if j==ny
                %vp1=0; %for lat=80
                %vn1=0;
                Fs=0; %no southward incoming at north face, from north
            else
                %vp1=(v(h,j+1,k)+abs(v(h,j+1,k)))/2; %for lat= -79 to 79
                vn1=(v(h,j+1,k)-abs(v(h,j+1,k)))/2;
                %UPWIND:
                %Fs=vn1*C(h,j+1,k);%southward incoming at north face, from box to north of j
                %CENTERED:
                Cf=(C(h,j,k)+C(h,j+1,k))/2;
                Fs=vn1*Cf*dxvall(j+1)*dzf(h);
            end
        
        % k loop: zonal velocity u- starts at west (0.5 centered lont(1))
            %along longitude- west and east faces- u and k
            %only need west face positive and east face negative vel.s:
            %%THIS IS PERIODIC! (for k=1, box 360 is possibly
            %%incoming)
            
            %TYPO: fixed. was up=(u(h,j,k)+abs(v(h,j,k)))/2;
            %at west face:
            if k==1
                %%eastward incoming at west face, from box to west of k
                up=(u(h,j,k)+abs(u(h,j,k)))/2; %1st k value 
                %UPWIND:
                %Fe=up*C(h,j,360); %eastward incoming at west face from box west of k
                %CENTERED:
                Cf=(C(h,j,360)+C(h,j,k))/2;
                Fe=up*Cf*dy*dzf(h);
            else
                up=(u(h,j,k)+abs(u(h,j,k)))/2; %eastward for lat = 0 to 359
                %UPWIND:
                %Fe=up*C(h,j,k-1); %eastward incoming at west face from box west of k
                %CENTERED:
                Cf=(C(h,j,k-1)+C(h,j,k))/2;
                Fe=up*Cf*dy*dzf(h);
            end
            
            %at east face:
            if k==nx
                %%westward flux from box to east (box 1) so k==1
                un1=(u(h,j,1)-abs(u(h,j,1)))/2; %negative vel means westward mvmt
                %UPWIND:
                %Fw=un1*C(h,j,1); %westward incoming at east face from the box east of k
                %CENTERED:
                Cf=(C(h,j,1)+C(h,j,k))/2;
                Fw=un1*Cf*dy*dzf(h);
            else
                un1=(u(h,j,k+1)-abs(u(h,j,k+1)))/2; %negative vel means westward mvmt
                %UPWIND:
                %Fw=un1*C(h,j,k+1); %westward incoming at east face from the box east of k
                %CENTERED:
                Cf=(C(h,j,k+1)+C(h,j,k))/2;
                Fw=un1*Cf*dy*dzf(h);
            end
            
        %h loop: vertical velocity w- starts at top (5m depth = zt(1))
            %from top to bottom, but here, POSITIVE w means UPWARD
            %velocity! also, i have 24 depths for w!!(but don't use 24)
            
            %OCCA update: w is same length as # of boxes, and starts at BOTTOM of box (w(1)=10m bottom of
            %box)
            
            if h==1
                Fd=0; %no flux from top face 
            else
                wn=(w(h-1,j,k)-abs(w(h-1,j,k)))/2; %top face,starting at 10m depth, NEGATIVE flux goes down into box 2
                %UPWIND:
                %Fd=wn*C(h-1,j,k);%downward flux from top face from box above h
                %CENTERED:
                Cf=(C(h-1,j,k)+C(h,j,k))/2;
                Fd=wn*Cf*dxall(j)*dy; %m3/s
            end
            
            %UPWARD flux is w POSITIVE
            if h==nz
                Fu=0;
            else
                 wp1=(w(h,j,k)+abs(w(h,j,k)))/2; %bottom face, starting with w(1) at 10m bottom of 1st box. POSITIVE is upward incoming
                 %UPWIND:
                 %Fu=wp1*C(h+1,j,k); %upward flux from box underneath h
                %CENTERED:
                Cf=(C(h+1,j,k)+C(h,j,k))/2;
                Fu=wp1*Cf*dxall(j)*dy;
            end
            
        %flux units: mmol/m3/s = m/s * mmol/m3 * 1/m
      
        %Feb 2017: for ref,
        %the distance between the faces (dzf, not dzc):
        %Fu=Fu*dxall(j)*dy; %depth of box 
        %Fd=Fd*dxall(j)*dy;
        %Fw=Fw*dy*dzf(h); %longitudinal distance of box (varies with lat)
        %Fe=Fe*dy*dzf(h);
        %Fn=Fn*dxvall(j)*dzf(h); %latitudinal distance of box
        %Fs=Fs*dxvall(j+1)*dzf(h);
        
        vol=dxall(j)*dy*dzf(h);
        
        advp(h,j,k)=(Fs-Fn+Fw-Fe+Fd-Fu)/vol; %positive fluxes are subtracted here.
        
        end
    end
end

end
