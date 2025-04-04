function hF2= bin_showFluorPerAtom(bindata,opts)
      
%% Get Data
% 
P=[bindata.Params];
X1 = [P.(opts.xVar)];

for kk=1:length(bindata)
    for rr=1:length(bindata(kk).LatticeBin)
        if isfield(bindata(kk).LatticeBin(rr),'PDF1_Center') && ~isempty(bindata(kk).LatticeBin(rr).PDF1_Center)

            Y(kk,rr)=bindata(kk).LatticeBin(rr).PDF1_Center;
            Ys(kk,rr)=bindata(kk).LatticeBin(rr).PDF1_Radius;
        else
            Y(kk,rr) = NaN;
            Ys(kk,rr) = NaN;
        end
    end
end


%% Special X-Variables

switch opts.xVar
    case 'qgm_Raman1_power'
        xVar2 = 'qgm_Raman2_power';
        X2 = [P.(xVar2)];
    case 'qgm_Raman2_power'
        xVar2 = 'qgm_Raman1_power';
        X2 = [P.(xVar2)];
    otherwise
        xVar2 = [];
        X2=[];
end
%% Plot it

hF2 = figure;
hF2.Color='w';
hF2.Position=[20 50 800 400];
hF2.Name = 'FluorPerAtom';


% ax =subplot(111);;  
ax=axes;
co=get(gca,'colororder');
legStr={};
markertype={'o','s','^','x'};

if ~isempty(xVar2)
    ux2 = unique(X2);
    myc = inferno(length(ux2));
    indsc=[];
    for nn=1:length(X2)
        indsc(nn) = find(X2(nn)==ux2,1);
    end
    cc = myc(indsc,:);
    for rr=length(bindata(1).LatticeBin):-1:1
        ps(rr)=scatter(X1,Y(:,rr),50,cc,markertype{rr},'filled');
        set(ps(rr),'MarkerEdgeColor','k')
    
        legStr{rr}=['image ' num2str(rr)];    hold on   
    end
    cb=colorbar;
    colormap(ax,inferno);
    caxis(ax,[min(X2) max(X2)]);
    cb.Label.String=xVar2;
    cb.Label.Interpreter='none';
else
    for rr=length(bindata(1).LatticeBin):-1:1
    ps(rr)=plot(X1,Y(:,rr),'o','linewidth',2,'markersize',6,'markerfacecolor',co(rr,:),...
        'markeredgecolor',co(rr,:)*.5);
    legStr{rr}=['image ' num2str(rr)];
    hold on   
    end
end
xlabel(opts.xVar,'interpreter','none')
if isequal(opts.xVar,'ExecutionDate')
    datetick x
end
ylabel('fluorescence/atom');
set(gca,'box','on','linewidth',1,'fontsize',10);
grid on;
yL=get(gca,'YLim');
set(gca,'YLim',[0 yL(2)]);

legend([ps],legStr,'location','best')


if isfield(bindata,'SourceDirectory') 
tFig=uicontrol('style','text','string',bindata(1).SourceDirectory,...
    'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','parent',hF2);
tFig.Position(4)=15;
tFig.Position(3)=hF2.Position(3);
tFig.Position(1:2)=[5 hF2.Position(4)-tFig.Position(4)];
end    


end

