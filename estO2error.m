%June 2020: EJZ 
%estimate error in oxygen concentrations
%ncdump('woa13_all_o00_01.nc')

o2_an = o_an;
o2_an(o2_an > 9e36)=nan;
o2_an = o2_an/22391.6*1000.0*1e3; %from ml/l to uM

o2_se = o_se; 
o2_se(o2_se > 9e36)=nan;
o2_se = o2_se/22391.6*1000.0*1e3; %from ml/l to uM
%nanmean(nanmean(nanmean(abs(o2_se)))) = 2.5 

%let's take the mean where O2 is lower than 50 uM:
o2_seL = o2_se;
o2_seL(o2_an>50)=nan;
%nanmean(nanmean(nanmean(abs(o2_seL)))) = 1.7 uM O2

%where O2 is lower than 10 uM:
o2_seL = o2_se;
o2_seL(o2_an>10)=nan;
%nanmean(nanmean(nanmean(abs(o2_seL)))) = 0.9 

%what about this standard error as a % uncertainty, on the grid?
uncO = o2_seL./o2_an;
uncOmean = nanmean(nanmean(nanmean(abs(uncO))));%= 0.12 %= 12% uncertainty below 10 uM
%Note: 0.06 uncertainty below 50 uM, and 0.02 overall, so uncertainty increases as O2 is
%lower -- use this 12% for the calculation

disp("Uncertainty of O2 due to error when concs are 10 uM or below is:")
disp(uncOmean)
disp("This is much larger than the uncertainty when O2 is high")

