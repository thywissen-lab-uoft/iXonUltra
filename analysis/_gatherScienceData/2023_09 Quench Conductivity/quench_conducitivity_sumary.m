d1=load('quench_1Er_2v2v_09_13');
d2=load('quench_1Er_4v4v_08_28');
d3 = load('quench_1Er_4v4v_09_12');

B2a = @(Bfield) 167*(1-(6.910)./(Bfield-202.15));


figure(2000);
clf
subplot(121);
co=get(gca,'colororder');
errorbar(d1.field,d1.N,d1.N_err,'o','markerfacecolor',co(1,:),'markersize',8,'linewidth',1,'markeredgecolor',co(1,:)*.5);
hold on
errorbar(d2.field,d2.N,d2.N_err,'o','markerfacecolor',co(2,:),'markersize',8,'linewidth',1,'markeredgecolor',co(2,:)*.5);
errorbar(d3.field,d3.N,d3.N_err,'o','markerfacecolor',co(3,:),'markersize',8,'linewidth',1,'markeredgecolor',co(3,:)*.5);
xlabel('field (G)');
ylabel('counts');


legend({'2v2v 9/13','4v4v 8/28','4v4v 9/12'});


% figure(2001);
% clf
% co=get(gca,'colororder');
subplot(122);
errorbar(B2a(d1.field).^2,d1.gamma,d1.gamma_err,'o','markerfacecolor',co(1,:),'markersize',8,'linewidth',1,'markeredgecolor',co(1,:)*.5);
hold on
errorbar(B2a(d2.field).^2,d2.gamma,d2.gamma_err,'o','markerfacecolor',co(2,:),'markersize',8,'linewidth',1,'markeredgecolor',co(2,:)*.5);
errorbar(B2a(d3.field).^2,d3.gamma,d3.gamma_err,'o','markerfacecolor',co(3,:),'markersize',8,'linewidth',1,'markeredgecolor',co(3,:)*.5);
xlabel('a^2 (a_0^2)');
ylabel('decay rate (1/ms)');
legend({'2v2v 9/13','4v4v 8/28','4v4v 9/12'});
