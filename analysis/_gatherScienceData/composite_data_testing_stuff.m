condata=[composite_data.conductivity_data(1)];
amp=[];
for kk=1:length(composite_data(1).conductivity_data)
   amp(kk)=[ composite_data(1).conductivity_data(kk).Params(1).conductivity_ODT2_mod_amp];
end



condata=[composite_data.conductivity_data];
Xs = [condata.XsBar];
XsErr = [condata.XsBarErr];
Ys = [condata.YsBar];
YsErr = [condata.YsBarErr];

C = [condata.C];
CErr=[condata.Cerr];
S = [condata.S];
SErr=[condata.Serr];

n=[condata.SpatialMaxUpDensity];
nerr=[condata.SpatialMaxUpDensityErr];
    

hF=figure(101);
hF.Color='w';
clf
co=get(gca,'colororder');

subplot(231);
errorbar(Xs,-C,CErr,CErr,XsErr,XsErr,'o','markerfacecolor',...
    co(1,:),'linewidth',1,'color',co(1,:)*.5);
xlabel('x sigma (um)');
ylabel('-C (um)');
ylim([-.5 2.5]);
title('\pi/2 phase amplitude');

subplot(232);
errorbar(Xs,S,SErr,SErr,XsErr,XsErr,'o','markerfacecolor',...
    co(2,:),'linewidth',1,'color',co(2,:)*.5);
xlabel('x sigma (um)');
ylabel('-S (um)');
ylim([-.5 2.5]);
title('0 phase amplitude');


subplot(233);
errorbar(Xs,n,nerr,nerr,XsErr,XsErr,'o','markerfacecolor',...
    co(3,:),'linewidth',1,'color',co(3,:)*.5);
xlabel('x sigma (um)');
ylabel('peak nup density');
ylim([0 0.15]);
title('density');

subplot(234);
total_amp = sqrt(C.^2+S.^2);
total_amp_err = (1./total_amp).*((C.*CErr).^2+(S.*SErr)).^(1/2);
errorbar(amp,total_amp,total_amp_err,'o','markerfacecolor',...
    co(4,:),'linewidth',1,'color',co(4,:)*.5);
xlabel('mod amp (V)');
ylabel('sqrt(C^2+S^2) (um)');
ylim([-.1 2.5]);
xlim([0 4.5]);
hold on
yyaxis right
amp_vec=linspace(0,4,3);
plot(amp_vec,amp_vec*2.63,'r-');
xlim([0 4.5]);
ylim([-.1 2.5]);
ylabel('drive amplitude (um)');
title('total amplitude and mod strength');

 subplot(235);
 errorbar(Xs,Ys,YsErr,YsErr,XsErr,XsErr,'o','markerfacecolor',...
    co(5,:),'linewidth',1,'color',co(5,:)*.5);
xlabel('x cloud size (um)');
ylabel('y cloud size (um)');
ylim([6 12]);
xlim([6 12]);
hold on
plot([3 12],[3 12],'k-');
axis equal tight
ylim([6 12]);
xlim([6 12]);
title('thermalization')


subplot(236);
errorbar((amp*2.63),total_amp./(amp*2.63),total_amp_err./(amp*2.63),'o','markerfacecolor',...
    co(4,:),'linewidth',1,'color',co(4,:)*.5);
xlabel('mod amp (um)');
ylabel('total amp/mod amp (um/um)');
ylim([0 2.5]);
 xlim([0 12]);
title('conductivity proxy')

