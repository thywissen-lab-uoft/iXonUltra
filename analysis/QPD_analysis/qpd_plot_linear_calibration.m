%% load data

qd = [ixondata.qpd_data];
P = [ixondata.Params];

X = [P.conductivity_ODT2_mod_amp];
X = [ixon_gaussdata.Xc(:,1)]'*16/(80);
Y = [qd.modfit_X1_C];

%% fit to line

 fit1=polyfit(X,Y,1);

%% plot

hF2 = figure;
hF2.Color='w';
hF2.Position=[1000 50 600 400];


    ax =axes;   
    
    co=get(gca,'colororder');

    hold on

  
    plot(X,Y,'o','linewidth',2,'markersize',10,'markerfacecolor',co(1,:),...
        'markeredgecolor',co(1,:)*.5);
    tVec = [min(X):0.01:max(X)];
    plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  
    
    
    xlabel('Gauss X center (um)','Interpreter','None')
    ylabel('ODT 1 QPD Norm. X (V/V)');
    set(gca,'box','on','linewidth',1,'fontsize',10);
    grid on;
    
legend('data','y='+string(fit1(1))+"x + "+string(fit1(2)))


%% load data

qd = [ixondata.qpd_data];
P = [ixondata.Params];

X = [P.conductivity_ODT2_mod_amp];
X = [ixon_gaussdata.Xc(:,1)]'*16/(80);
Y = [qd.modfit_X2_C];

%% fit to line

 fit1=polyfit(X,Y,1);

%% plot

hF2 = figure;
hF2.Color='w';
hF2.Position=[1000 50 600 400];


    ax =axes;   
    
    co=get(gca,'colororder');

    hold on

  
    plot(X,Y,'o','linewidth',2,'markersize',10,'markerfacecolor',co(1,:),...
        'markeredgecolor',co(1,:)*.5);
    tVec = [min(X):0.01:max(X)];
    plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  
    
    
    xlabel('Gauss X center (um)','Interpreter','None')
    ylabel('ODT 2 QPD Norm. X (V/V)');
    set(gca,'box','on','linewidth',1,'fontsize',10);
    grid on;
    
legend('data','y='+string(fit1(1))+"x + "+string(fit1(2)))
      