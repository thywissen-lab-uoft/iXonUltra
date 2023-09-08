function hF=ixon_showSize(data,xVar,plt_opts,fit_opts)

if nargin <4 
    fit_opts = struct;
end

if nargin < 3 
    plt_opts = struct;
end

if nargin < 2
    xVar = 'ExecutionDate';
end

if ~isfield(plt_opts,'FigLabel')
    plt_opts.FigLabel = '';
end
%% Grab Data

params=[data.Params];
xvals=[params.(xVar)];

PixelSize = data.PixelSize;
Xs = data.Xs*PixelSize/data.Magnification;
Ys = data.Ys*PixelSize/data.Magnification;

%% Make Figure

hF=figure('Name',[pad(['ixon' data.FitType ' size'],20) plt_opts.FigLabel],...
    'units','pixels','color','w','numbertitle','off');
hF.Position(1)=1015;
hF.Position(2)=50;
hF.Position(3)=800;
hF.Position(4)=300;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',plt_opts.FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','iXon','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax1=subplot(131);
set(hax1,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' plt_opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(Xs,2)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,Xs(:,nn),'o-','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end


str=[data.FitType newline ' $\sigma_X (\mu \mathrm{m})$'];
text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


% Make axis

hax2=subplot(132);
set(hax2,'box','on','linewidth',1,'fontsize',10,...
    'ygrid','on','xgrid','on');
hold on
xlabel([xVar ' (' plt_opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(Ys,2)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,Ys(:,nn),'o-','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end


str=[data.FitType newline ' $\sigma_Y (\mu \mathrm{m})$'];
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


% Make axis
hax3=subplot(133);
set(hax3,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' plt_opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(Xs,2)
    [cface,cedge] = ixoncolororder(nn);
   plot(xvals,pi*Xs(:,nn).*Ys(:,nn),'o-','color',cedge,'linewidth',1,'markersize',8,...
       'markerfacecolor',cface,'markeredgecolor',cedge);
end

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

str=[data.FitType newline '$\pi \sigma_X \sigma_Y (\mu \mathrm{m}^2)$'];
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');

% hax1.Position(4)=hax1.Position(4)-15;
% hax2.Position(4)=hax1.Position(4);
% hax3.Position(4)=hax1.Position(4);

ixon_resizeFig(hF,t,[hax1 hax2 hax3]);


