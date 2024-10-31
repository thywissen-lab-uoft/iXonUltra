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
hF.Position= [100 100 1200 600];
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

subplot(231)
plot(digdata.X,digdata.Xc_um,'o','markerfacecolor',co(1,:),...
'linewidth',1,'markeredgecolor',co(1,:)*.5);
ylabel('x position (\mum)');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end


subplot(234)
plot(digdata.X,digdata.Xs_um,'o','markerfacecolor',co(1,:),...
'linewidth',1,'markeredgecolor',co(1,:)*.5);
ylabel('x sigma (\mum)');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end

subplot(232)
plot(digdata.X,digdata.Yc_um,'o','markerfacecolor',co(2,:),...
'linewidth',1,'markeredgecolor',co(2,:)*.5);
ylabel('y position (\mum)');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end

subplot(235)
plot(digdata.X,digdata.Ys_um,'o','markerfacecolor',co(2,:),...
'linewidth',1,'markeredgecolor',co(2,:)*.5);
ylabel('y sigma (\mum)');
xlabel(digdata.xVar,'interpreter','none');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end

subplot(233)
plot([digdata.X],[digdata.Natoms],'ko','markerfacecolor',[.5 .5 .5],...
'linewidth',1);
xlabel(digdata.xVar,'interpreter','none');
ylabel('Natoms');
if isequal(digdata.xVar,'ExecutionDate')
    datetick x
end
if isfield(digdata,'ThresholdingType')
    if digdata.ThresholdingType == 'CompensatedInd'
        str = [digdata.ThresholdingType ' thresholding' newline '$N_{\mathrm{thresh}} = ' num2str(mean([digdata.Threshold])) '$'];
    else  
        str = [digdata.ThresholdingType ' thresholding' newline '$N_{\mathrm{thresh}} = ' num2str(unique([digdata.Threshold])) '$'];
    end
else
    str = ['$N_{\mathrm{thresh}} = ' num2str(unique([digdata.Threshold])) '$'];
end
text(2,2,str,'units','pixels','fontsize',10,'verticalalignment','bottom',...
'horizontalalignment','left','interpreter','latex');

subplot(236)

Nbar = mean([rmoutliers(digdata.Natoms)]);
Nstd = std([rmoutliers(digdata.Natoms)]);

histogram(rmoutliers(digdata.Natoms),15);
  xlabel('atom number');
ylabel('occurences');
text(1,1,'*outliers removed*','units','normalized','fontsize',8,'verticalalignment','cap',...
'horizontalalignment','right');
str = ['N = ' num2str(round(Nbar)) ' \pm ' num2str(round(Nstd))];
text(.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
'horizontalalignment','left','backgroundcolor',[1 1 1 .5]);


end

