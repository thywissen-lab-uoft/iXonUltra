function [outputArg1,outputArg2] = dig_showSummaryAverage(digdata,opts)

if nargin ==1
    opts=struct;
end

if ~isfield(opts,'doSmooth')
    opts.doSmooth=false;
end

if ~isfield(opts,'smoothRadius')
    opts.SmoothRadius = 5;
end

%% Initialize Figure

if ~isfield(opts,'Parent')
    opts.Parent = figure('color','w','Position',[100 100 800 600],...
        'Name','DigSummary','NumberTitle','off');
    fig = opts.Parent;
else
    fig = opts.Parent;
    for kk=1:length(fig.Children)
        delete(fig.Children(1))
    end
end
    
if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=400;
    tFig.Position(1:2)=[1 1];
end    

%% Average Image Calculation

z_bar = mean(digdata.Zdig(:,:,1),3);
nImages = size(digdata.Zdig,3);

z_1    = sum(z_bar,1)/size(digdata.Zdig,3);z_1=z_1(:);
z_2    = sum(z_bar,2)/size(digdata.Zdig,3);z_2=z_2(:);

n1 = digdata.n1;
n2 = digdata.n2;

radial_step = 3;

 [Tics,Average,dev,n]=radial_profile(z_bar,radial_step)

%% Main Image

strDesc = [num2str(nImages) ' images'];

axDig = subplot(1,3,1,'parent',fig);
imagesc(n1,n2,z_bar,'parent',axDig);
axis(axDig,'equal');
axis(axDig,'tight');
set(axDig,'box','on','linewidth',1,'ydir','normal','xaxislocation','top',...
    'yaxislocation','left','fontsize',8);
xlabel(axDig,'site 1');
ylabel(axDig,'site 2');
colormap(axDig,'bone');
% cc=colorbar(axDig,'location','west');
% cc.Label.String='occupation';
% cc.Label.FontSize=6;
clim(axDig,[0 .2]);

if nImages>1
    text(.01,.01,strDesc,'fontsize',6,'color','r','HorizontalAlignment','left',...
        'verticalalignment','bottom','units','normalized');
end

%% n2 plot
% ax2 = subplot(2,3,2,'parent',fig);
% plot(z_2,n2,'-')
% 
% %% n1 plot
% ax1 = subplot(2,3,4,'parent',fig);
% plot(n1,z_1,'-')

%% radial plot
axr = subplot(1,3,2,'parent',fig);
plot(Tics,Average)
end


function [Tics,Average,dev,n]=radial_profile(data,radial_step)
%main axii cpecified:
x=(1:size(data,2))-size(data,2)/2;
y=(1:size(data,1))-size(data,1)/2;
% coordinate grid:
[X,Y]=meshgrid(x,y);
% creating circular layers
Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % illustrating the principle:
% % figure;imagesc(Z_integer.*data)
% very fast MatLab calculations:
Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
Average=accumarray(Z_integer(:),data(:),[],@mean);


% dev=accumarray(Z_integer(:),data(:),[],@std);
dev=accumarray(Z_integer(:),data(:),[],@(x) std(x,1));

n= accumarray(Z_integer(:),data(:),[],@(x) numel(x));

end


function [axX,axY,axWidth,axHeight]=getAxesPos(num_rows,num_cols,index,figW,figH)
yTop=30;
yBot=50;

xLeft=35;
xRight=20;

ySpace=35;
xSpace=70;

rowNumber = floor((index-1)/num_cols)+1;
colNumber = mod(index-1,num_cols)+1;

axHeight = (figH-yTop-yBot-ySpace*(num_rows-1))/num_rows;
axWidth = (figW-xLeft-xRight-xSpace*(num_cols-1))/num_cols;

axX = xLeft + (axWidth+xSpace)*(colNumber-1);
axY=(figH-yTop-axHeight)-(rowNumber-1)*(axHeight+ySpace);

end

