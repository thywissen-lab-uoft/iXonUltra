function hF = dig_showFidelitySummary(digdata,opts)

if nargin~=2
    opts=struct;
end
%% Analysis

pdLoss = fitdist([digdata.lost_fraction]','normal');
Rloss_mean=pdLoss.mu; 
Rloss_sigma=pdLoss.sigma;
loss_str = ['$' num2str(round(100*Rloss_mean,1)) '\% \pm' num2str(round(100*Rloss_sigma,1)) '\% $'];

pdLossCenter = fitdist([digdata.lost_fraction_center]','normal');
RlossCenter_mean=pdLossCenter.mu; 
RlossCenter_sigma=pdLossCenter.sigma;
losscen_str = ['$' num2str(round(100*RlossCenter_mean,1)) '\% \pm' num2str(round(100*RlossCenter_sigma,1)) '\% $'];

pdHop = fitdist([digdata.hop_fraction]','normal');
Rhop_mean=pdHop.mu; 
Rhop_sigma=pdHop.sigma;
hop_str = ['$' num2str(round(100*Rhop_mean,1)) '\% \pm' num2str(round(100*Rhop_sigma,1)) '\% $'];

% Peak charge density if it were a gaussian
nPeakGauss=mean(digdata.nPeakGauss,2);

%% Make the figure
if ~isfield(opts,'Parent') || isempty(opts.Parent)
    opts.Parent = figure;
    set(opts.Parent,'color','w','Name',['fidelity_summary']);
    clf

    opts.Parent.Position=[50 50 700 500];
end

hF=opts.Parent;

if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
    t=uicontrol('style','text','string',opts.FigLabel,'fontsize',7,...
        'backgroundcolor','w','Position',[1 1 600 15],'horizontalalignment','left',...
        'parent',hF);
    % t.Position(2) = t.Parent.Position(4)-t.Position(4)-2;
end
 
ax1 = subplot(211);
co=get(gca,'colororder');
pL_all=plot([digdata.X],100*[digdata.lost_fraction],'o','markerfacecolor',co(1,:),...
    'color',co(1,:)*.5,'linewidth',1);
hold on
pL_cen=plot([digdata.X],100*[digdata.lost_fraction_center],'o','markerfacecolor',co(2,:),...
    'color',co(2,:)*.5,'linewidth',1);
 
pH=plot([digdata.X],100*[digdata.hop_fraction],'^','markerfacecolor',co(3,:),...
    'color',co(3,:)*.5,'linewidth',1);


xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x;
end
ylim([0 30]);
ylabel('rate (%)');

yyaxis right
set(gca,'YColor',[.5 .5 .5]);
pDensity=plot([digdata.X],nPeakGauss,'x','color',[.5 .5 .5]);
ylim([0 .5])
ylabel('charge density')

legend([pL_all,pL_cen,pH,pDensity],...
    {'loss (all R) ',...
    ['loss (R<' num2str(digdata.FidelityCenterRadius) ')'],...
    'hop (all R)','peak gauss'},'orientation','horizontal','location','northwest')

subplot(234);
p1=histfit([digdata.lost_fraction]);
title('total loss fraction')
ylabel('occurences')
xlabel('rate');
legend(p1(2),loss_str,'interpreter','latex')

subplot(235);
p2=histfit([digdata.lost_fraction_center]);
title(['center loss fraction (R<' num2str(digdata.FidelityCenterRadius) ')']);
ylabel('occurences')
xlabel('rate');
legend(p2(2),losscen_str,'interpreter','latex')

subplot(236);
p3=histfit([digdata.hop_fraction]);
title('hop fraction')
ylabel('occurences')
xlabel('rate');
legend(p3(2),hop_str,'interpreter','latex')

end
