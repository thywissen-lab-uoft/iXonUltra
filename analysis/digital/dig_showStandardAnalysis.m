function hF =  dig_showStandardAnalysis(digdata,opts)



if nargin ==1
    opts = struct;
end

if ~isfield(opts,'doAnimate')
    opts.doAnimate = 0;
end

if ~isfield(opts,'doSave')
    opts.doSave = 0;
end

if ~isfield(opts,'ROI')
    opts.ROI = 'max';
end

%% Get ROI
    
  

%% Initialize Graphics
hF = figure;
clf
hF.Color='w';
hF.Position= [100 100 1200 500];
hF.Name = 'Digital Standard';

if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=hF.Position(3);
    tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
end    
co=get(gca,'colororder');

subplot(241)
plot(digdata.X,digdata.Xc_um,'o','markerfacecolor',co(1,:),...
'linewidth',1,'markeredgecolor',co(1,:)*.5);
ylabel('mean(X) (\mum)');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end


subplot(245)
plot(digdata.X,digdata.Xs_um,'o','markerfacecolor',co(1,:),...
'linewidth',1,'markeredgecolor',co(1,:)*.5);
ylabel('\sigma_x (\mum)');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end

subplot(242)
plot(digdata.X,digdata.Yc_um,'o','markerfacecolor',co(2,:),...
'linewidth',1,'markeredgecolor',co(2,:)*.5);
ylabel('mean(Y) (\mum)');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end

subplot(246)
plot(digdata.X,digdata.Ys_um,'o','markerfacecolor',co(2,:),...
'linewidth',1,'markeredgecolor',co(2,:)*.5);
ylabel('\sigma_y (\mum)');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end

subplot(243)
plot([digdata.X],[digdata.Natoms],'ko','markerfacecolor',[.5 .5 .5],...
'linewidth',1);
xlabel(digdata.xVar,'interpreter','none');
ylabel('charge');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end
% Nbar = mean([rmoutliers(digdata.Natoms)]);
% Nstd = std([rmoutliers(digdata.Natoms)]);
Nbar = mean(digdata.Natoms);
Nstd = std(digdata.Natoms);
Nmed = median(digdata.Natoms);

Nstr = ['$\mathrm{mean}(N) : ' num2str(round(Nbar)) '$' newline ...
    '$\mathrm{med}(N) : ' num2str(round(Nmed)) '$' newline ...
    '$\mathrm{std}(N) : ' num2str(round(Nstd)) '$'];
text(.01,.99,Nstr,'units','normalized','fontsize',7,...
    'interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left');


if isfield(digdata,'ThresholdingType')
    if isequal(digdata.ThresholdingType,'CompensatedInd')
        str = ['compensated threshold'];
    end     
else
    str = ['$\mathrm{thresh} = ' num2str(unique([digdata.Threshold])) '$'];
end

text(0.01,.01,str,'units','normalized','fontsize',8,...
    'horizontalalignment','left','verticalalignment','bottom')

subplot(247)
plot([digdata.X],[digdata.npeak],'o','markerfacecolor',co(6,:),'markeredgecolor',0.5*co(6,:),...
        'linewidth',2,'color',co(6,:)*.5);
ylabel('peak gauss density n_{0,\uparrow}');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end
xlabel(digdata.xVar,'interpreter','none');
str = '$n_{0,\uparrow}=0.5N/(2\pi\sigma_x\sigma_y)$';
text(.01,.99,str,'units','normalized','fontsize',10,'verticalalignment','top',...
'horizontalalignment','left','interpreter','latex');
% 
% 
% subplot(236)
% 
% Nbar = mean([rmoutliers(digdata.Natoms)]);
% Nstd = std([rmoutliers(digdata.Natoms)]);
% 
% histogram(rmoutliers(digdata.Natoms),15);
%   xlabel('atom number');
% ylabel('occurences');
% text(1,1,'*outliers removed*','units','normalized','fontsize',8,'verticalalignment','cap',...
% 'horizontalalignment','right');
% str = ['N = ' num2str(round(Nbar)) ' \pm ' num2str(round(Nstd))];
% text(.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
% 'horizontalalignment','left','backgroundcolor',[1 1 1 .5]);


end

