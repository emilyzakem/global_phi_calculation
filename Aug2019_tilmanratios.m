%FIGURE: SUPPLY in RR space
figure;
tfont=16;

Dinline=0:0.01:.8; %mock POCdiv
philine=Dinline*r;

logo2=o2conc;
logo2(logo2<=0)=nan;
logo2=log10(logo2);

E = philine*.35; %error from POC flux

plot(Dinline,philine,'k',Dinline,[philine+E; philine-E],'k--','linewidth',2)
%this one shows the red lines and the blue lines: y and POC error:
%plot(Dinline,philine,'k',Dinline,philine+E,'b',Dinline,philine-E,'b--',Dinline,Dinline*rU,'r',Dinline,Dinline*rL,'r--','linewidth',2)


hold on;
scatter(POCdiv(:),O2in(:),10,logo2(:),'filled') %M
%scatter(POCdiv(:),O2in(:),10,o2conc(:),'filled') %M
%scatter(Dinp(:),O2inp(:),25,pbnp(:),'filled')

%just black lines:
plot(Dinline,philine,'k',Dinline,[philine+E; philine-E],'k--','linewidth',2)
%with red and blue:
%plot(Dinline,philine,'k',Dinline,philine+E,'b',Dinline,philine-E,'b--',Dinline,Dinline*rU,'r',Dinline,Dinline*rL,'r--','linewidth',2)

%LOG axes:
set(gca,'yscale','log')
%set(gca,'xscale','log')
ylabel('O_2_i_n (\muM day^-^1) ','fontsize',16)
xlabel('OM_d_i_v (\muM POC day^-^1) ','fontsize',16)
title('O_2 (log M)','fontsize',16)
%just black lines:
legend('\phi = 1 ', 'OM_i_n 35% error' )%,'','Open ocean lower y')
%with red and blue:
%legend('\phi = 1 ', 'OM_i_n+35%','OM_i_n-35%','r for BGE=0.35','r for BGE=0.02' )%,'','Open ocean lower y')
colorbar

caxis([0 2.7])

%xlim([-0.02 .3])
xlim([0 .3])
%ylim([1e-3 1.5e2])
ylim([1e-3 300])

set(gcf,'color','w');
set(gca,'fontsize',tfont);