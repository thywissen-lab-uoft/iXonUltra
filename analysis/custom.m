% Center frequency for expected RF field (if relevant)
B = atomdata(1).Params.HF_FeshValue_Initial;
x0= (BreitRabiK(B,9/2,-5/2)-BreitRabiK(B,9/2,-7/2))/6.6260755e-34/1E6;


% Grab Raw data
X=Ndatabox.X;   
X=X-x0;    
X=X*1E3;  
X=X';
xstr=['frequency - ' num2str(round(abs(x0),4))  ' MHz (kHz)'];

% Define Y Data
N1=Ndatabox.Natoms(:,1);
N2=Ndatabox.Natoms(:,2);     
N2=N2/0.6;

Y=N1+N2;


hF=figure;
hF.Position(3:4)=[800 300];
set(hF,'color','w','name','custom box');
co=get(gca,'colororder');

subplot(131);
plot(X,N1,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
    'linewidth',2,'markersize',8);
xlabel('$\Delta f~(\mathrm{kHz})$','interpreter','latex');
ylabel('N_{9/2}');    
ylim([0 1.2E5]);
set(gca,'fontsize',12,'linewidth',1,'box','on','xgrid','on','ygrid','on');


subplot(132);
plot(X,N2,'o','markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5,...
    'linewidth',2,'markersize',8);
xlabel('$\Delta f~(\mathrm{kHz})$','interpreter','latex');
ylabel('N_{7/2}');    
ylim([0 1.2E5]);
set(gca,'fontsize',12,'linewidth',1,'box','on','xgrid','on','ygrid','on');

subplot(133);
plot(X,N1+N2,'o','markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.5,...
    'linewidth',2,'markersize',8);
xlabel('$\Delta f~(\mathrm{kHz})$','interpreter','latex');
ylabel('N_{tot}');    
ylim([0 2.2E5]);
set(gca,'fontsize',12,'linewidth',1,'box','on','xgrid','on','ygrid','on');

