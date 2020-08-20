%Modified (annotated) from O2andPOC_Top_Aug2019.m by EJZ
%Jun 14, 2020
%Further improvements -- got all 4 components of diff working on Jun 24, 2020


%load in O2in and POCdiv
disp('O2in loading...')
O2in_Jun2020_Load
    
disp('O2in loaded. POCdiv loading...')
POCdiv_Load

disp('POCdiv loaded. O2 interpolating...')
O2in_Interp
          
disp('done')
       
%now have:
%O2in in uM O2/day
%POCin in uM C/day
%NO3in in uM NO3/day (check this)

o2conc=O2interp;% in uM O2 %concentration
lat=latosh;
lon=lono;
zc; %m, centered depth
dx; %m, 160 long 
dy; %m, =111000
dzf; %dz at faces


%%
%Plot Phi:

r=1/.9

rObs=O2in./POCdiv;

phi = rObs/r;

%new uncertainty, added June 2 2020: 
%from https://en.wikipedia.org/wiki/Propagation_of_uncertainty, assuming NO
%covariance (there actually IS some covariance in O2 and POC but that's not
%quantifiable, i.e. i don't know what it is or how to calc it)
%stdf = abs(f)*sqrt( (std1/f1)^2 + (std2/f2)^2 + (std3/f3)^3 )

%for f +- err

%errR = 0.1*r; %10%error
%errO = 0.12*O2in; %X%error O2 from estO2error_Jun2020.m
%halfway betw 20 and 35% for POC flux: 27.5
uncert = sqrt( 0.275^2 + 0.1^2 + 0.12^2 ); % = 0.32

err = phi*uncert;
disp('Propogated uncertainty is: ')
disp(uncert)

phiU = phi + err;
phiL = phi - err;

%%
%send files to .mat for plotting in python:


%saved Jun 14 2020 with 27.5% uncertainty in POC flux, 32% uncert overall:
%save('/Users/ezakem/Documents/DeepRespiration/doc/Paper_DiagnosticPhi/Figures/VarsforMapPlot_061420.mat',...
%    'lat','lon',...
%    'O2in','POCdiv','o2conc','ONESin','phi','zc','O2','zo','lato','lono',...
%    'uncert','phiU','phiL')

%saved updated Jun 24 2020, with implicit vert diff and ALL 4 components of
%mixing now correctly calculated to give positive incoming fluxes only (min
%0 for all)
%save('/Users/ezakem/Documents/DeepRespiration/doc/Paper_DiagnosticPhi/Figures/VarsforMapPlot_062420.mat',...
%    'lat','lon',...
%    'O2in','POCdiv','o2conc','ONESin','phi','zc','O2','zo','lato','lono',...
%    'uncert','phiU','phiL') %Nloss = 42.1

%saved updated Jun 25 2020, with UPDATED O2in_diffp calcs in toGithub_v2 with new
%save('/Users/ezakem/Documents/DeepRespiration/doc/Paper_DiagnosticPhi/Figures/VarsforMapPlot_062520.mat',...
%    'lat','lon',...
%    'O2in','POCdiv','o2conc','ONESin','phi','zc','O2','zo','lato','lono',...
%    'uncert','phiU','phiL') %Nloss = 42.3
%EXPLICIT z2 (gives similar values as implicit, but more transparent -- see
%email to Jon)
%save('/Users/ezakem/Documents/DeepRespiration/doc/Paper_DiagnosticPhi/Figures/VarsforMapPlot_062520_explz2.mat',...
%    'lat','lon',...
%    'O2in','POCdiv','o2conc','ONESin','phi','zc','O2','zo','lato','lono',...
%    'uncert','phiU','phiL') %Nloss = 42.5

%%

plotd=4; %with ref to zc

figure;
phid=squeeze(phi(plotd,:,:));

o2d=squeeze(o2conc(plotd,:,:));
%ONESdivd=squeeze(ONESdiv(plotd,:,:));
%onessumrd=squeeze(onessumr(plotd,:,:));
%get abs value, since O2 div can be pos and neg.
%rObsd=abs(abs(rObsd));
%or, take out negative:
%rObsd(rObsd<0)=2; %blue on fringe..
%and/or, take out where o2conc > X uM:
%rObsd(o2d>100)=2;
pcolor(lon,lat,phid)
%hold on; contour(lon,lat,rObsd,[1 1],'k','linewidth',2); hold off

%pcolor(lon,lat,ONESdivd) %there is some divergence
%pcolor(lon,lat,onessumrd*100) 
shading flat;colorbar
title(['\Phi (mol O2/mol C*r) at ' num2str(zc(plotd)) 'm'])
%title(['\Phi, removing \Phi<0, masking O2>100 \muM '])

%caxis([0 300])
caxis([0 2])
xlim([0 360])
ylim([-80 80])
makerbcolormap;
set(gcf,'color','w');
set(gca,'fontsize',14);
%xlabel('Longitude ')
ylabel('Latitude ')
set(gca,'xtick',[])

%%

Aug2019_tilmanratios

%%

%turn the following off and on to get uncertainty in Nloss
%phi = phiU;
%phi = phiL;

Calc_Nloss




