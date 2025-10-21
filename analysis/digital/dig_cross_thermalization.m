function [hF_XT,dig_XT_data] = dig_cross_thermalization(digdata,opts)
%% Load digdata
d = digdata;


%% Recalculate 2nd moment in ODT basis

a_px = digdata.Lattice_px;
a_um = digdata.Lattice_um;  

for nn=1:length(digdata.FileNames)
        for rr=1:size(digdata.Zdig,4)
            Ratom = digdata.Ratom{nn,rr};
            
            Xs_px   = 0.5*std(Ratom(1,:)*sqrt(3)-Ratom(2,:));
            Xs_site = Xs_px/a_px;
            Xs_um(nn)   = Xs_site*a_um;
        
            Ys_px   = 0.5*std(Ratom(1,:)+sqrt(3)*Ratom(2,:));
            Ys_site = Ys_px/a_px;
            Ys_um(nn)   = Ys_site*a_um;       
            
            nPeakGauss(nn,rr) = digdata.Natoms(nn,rr)./(sqrt(2*pi*Xs_site.^2).*sqrt(2*pi*Ys_site.^2));
              
        end  
end

X = d.X;

Xs = d.Xs_um;
Ys = d.Ys_um;

Xs = Xs_um';
Ys = Ys_um';

ratio = Xs./Ys;
% Y = Ys;

% keyboard

%% Exponential fit

myfit=fittype('A*exp(-t/tau)+offset','coefficients',{'A','tau','offset'},...
    'independent','t');

%% Plot and fit the data

hF_XT = figure(1000);
clf;
hF_XT.Color='w';
hF_XT.Position=[10 50 800 600];
hF_XT.Name = 'HeatingRate';
co=get(gca,'colororder');

if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=hF_XT.Position(3);
    tFig.Position(1:2)=[5 hF_XT.Position(4)-tFig.Position(4)];
end    
%%%%%  Xs %%%%
axX = subplot(2,2,1);

% Plot data
plot(X,Xs,'ko','markersize',8,'markerfacecolor','r','LineWidth',1)
hold on

% Fit 
opt=fitoptions(myfit);
opt.StartPoint = [Xs(1)-Ys(end),max(X)/2, Xs(end)]; 
opt.Lower = [-Inf -Inf 1];
foutX=fit(X',Xs,myfit,opt);

Xplot = [0:1:max(X)];
plot(Xplot,feval(foutX,Xplot),'k-','linewidth',2)

xlabel('Cross-thermalization hold time (ms)')
ylabel('\sigma_x (um)')
str = ['y = ' num2str(foutX.A,2) '*exp(-t/' num2str(foutX.tau,3) 'ms) + ' num2str(foutX.offset,3) ];
legend('Data',str,'Location','best')

%%%%%  Ys %%%%
axY = subplot(2,2,2);

% Plot data
plot(X,Ys,'ko','markersize',8,'markerfacecolor','r','LineWidth',1)
hold on

% Fit 
opt=fitoptions(myfit);
opt.StartPoint = [Ys(1)-Ys(end),max(X)/2, Ys(end)];   
opt.Lower = [-Inf -Inf 1];
foutY=fit(X',Ys,myfit,opt);

Xplot = [0:1:max(X)];
plot(Xplot,feval(foutY,Xplot),'k-','linewidth',2)

xlabel('Cross-thermalization hold time (ms)')
ylabel('\sigma_y (um)')
str = ['y = ' num2str(foutY.A,2) '*exp(-t/' num2str(foutY.tau,3) 'ms) + ' num2str(foutY.offset,3) ];
legend('Data',str,'Location','best')

%%%%% Ratio %%%%
axR = subplot(2,2,3);

% Plot data
plot(X,ratio,'ko','markersize',8,'markerfacecolor','r','LineWidth',1)
hold on

% Fit 
opt=fitoptions(myfit);
opt.StartPoint = [ratio(1)-1,max(X)/2, ratio(end)];  
opt.Lower = [-Inf -Inf 0.5];
opt.Upper = [Inf Inf 1.5];
foutR=fit(X',ratio,myfit,opt);

Xplot = [0:1:max(X)];
plot(Xplot,feval(foutR,Xplot),'k-','linewidth',2)

xlabel('Cross-thermalization hold time (ms)')
ylabel('\sigma_x/\sigma_y')
str = ['y = ' num2str(foutR.A,2) '*exp(-t/' num2str(foutR.tau,3) 'ms) + ' num2str(foutR.offset,3) ];
legend('Data',str,'Location','best')

%%%%% Energy %%%%
axR = subplot(2,2,4);

% Plot data
plot(X,Xs.^2+Ys.^2,'ko','markersize',8,'markerfacecolor','r','LineWidth',1)
hold on

% % Fit 
% opt=fitoptions(myfit);
% opt.StartPoint = [ratio(1)-1,max(X)/2, ratio(end)];  
% opt.Lower = [-Inf -Inf 1];
% opt.Upper = [Inf Inf 1];
% foutR=fit(X',ratio,myfit,opt);
% 
% Xplot = [0:1:max(X)];
% plot(Xplot,feval(foutR,Xplot),'k-','linewidth',2)

xlabel('Cross-thermalization hold time (ms)')
ylabel('\sigma_x^2 + \sigma_y^2 (um^2)')
% str = ['y = ' num2str(foutR.A,2) '*exp(-t/' num2str(foutR.tau,3) 'ms) + ' num2str(foutR.offset,3) ];
% legend('Data',str,'Location','best')

%% Output XT_data
clear dig_XT_data
dig_XT_data            = struct;
dig_XT_data.fitX       = foutX;
dig_XT_data.fitY       = foutY;
dig_XT_data.fitR       = foutR;
dig_XT_data.ratio      = ratio;
dig_XT_data.nPeakGauss = nPeakGauss;















end

