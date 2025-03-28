function hF2= bin_showFluorPerAtom(bindata,opts)
      
%% Get Data
% 
P=[bindata.Params];
X = [P.(opts.xVar)];

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
%% Plot it

hF2 = figure;
hF2.Color='w';
hF2.Position=[1000 50 600 400];
hF2.Name = 'FluorPerAtom';
ax =axes;       
co=get(gca,'colororder');
legStr={};
    
for rr=1:length(bindata(1).LatticeBin) 
    ps(rr)=errorbar(X,Y(:,rr),2*Ys(:,rr),'o','linewidth',2,'markersize',6,'markerfacecolor',co(rr,:),...
        'markeredgecolor',co(rr,:)*.5);
    legStr{rr}=['image ' num2str(rr)];
    hold on   
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

legend([ps],legStr)


if isfield(bindata,'SourceDirectory') 
tFig=uicontrol('style','text','string',bindata(1).SourceDirectory,...
    'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','parent',hF2);
tFig.Position(4)=15;
tFig.Position(3)=hF2.Position(3);
tFig.Position(1:2)=[5 hF2.Position(4)-tFig.Position(4)];
end    


end

