function bin_Fidelity(data,opts)

if nargin ==1
    opts =struct;
end

if ~isfield(opts,'threshold')
    opts.threshold=1500;
end

if ~isfield(opts,'FigureNumber')
   opts.FigureNumber=20041; 
end

Nthresh = opts.threshold;           

 ca = [0 0 0];       
cb = [0.7 .1 .6];
cc = [linspace(ca(1),cb(1),1000)' ...
    linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
                
           %%     
hF = figure(opts.FigureNumber);
set(hF,'color','w','Name','Bin Fidelity');
clf
hF.Position=[800 100 1100 800];

if isfield(opts,'Label') && ~isempty(opts.Label)
uicontrol('style','text','string',opts.Label,'fontsize',7,...
    'backgroundcolor','w','Position',[1 1 hF.Position(3) 15],'horizontalalignment','center');
end

Z1 = data.LatticeBin(1).Zbin;
Z1b = Z1;
Z1b(isnan(Z1)) = 0;
N1 = sum(Z1b,'all');

subplot(2,2,1);
imagesc(data.LatticeBin(1).n1,data.LatticeBin(1).n2,Z1);
colormap(cc);
axis equal tight
title('image 1');
caxis([0 2*Nthresh]);
title('image 1');
set(gca,'ydir','normal','box','on','linewidth',1);

ax_hB1=subplot(2,2,2);
pHistB1 = bar(1:100,1:100,'parent',ax_hB1,'linestyle','none',...
    'facecolor','k');
hold on
pHistBdivide = plot([1 1]*50,[0 100],'k-','parent',ax_hB1);
ylabel('occurences');
xlabel('counts per lattice site');
set(ax_hB1,'box','on','linewidth',.1,'fontsize',12,'units','normalized',...
    'XAxisLocation','bottom','YDir','normal','UserData','H1');
yyaxis right
pHistB2 = bar(1:100,1:100,'parent',ax_hB1,'linestyle','none',...
    'FaceColor',[0.6 0 0.5]);
ylabel('occurences');
ax_hB1.YColor=[0.6 0 0.5];   
x = data.LatticeHistogram(1).Centers;
xe = data.LatticeHistogram(1).Edges;
y = data.LatticeHistogram(1).N; 
xL = x<=Nthresh;
xH = ~xL;        
set(pHistB1,'XData',x(xL),'YData',y(xL));
set(pHistB2,'XData',x(xH),'YData',y(xH));        
pHistBdivide.Parent.YAxis(1).Limits = [0 max(pHistB1.YData)*1.1];
pHistBdivide.Parent.YAxis(2).Limits = [0 max(pHistB2.YData)*1.1];        
pHistBdivide.Parent.XAxis.Limits = [0 max(x)*1.1];
set(pHistBdivide,'Xdata',[1 1]*Nthresh,'Ydata',pHistBdivide.Parent.YAxis(1).Limits);
title('image 1');
Nlow = sum(y(xL));
Nhigh = sum(y(xH));
str = ['threshold : ' num2str(Nthresh) newline ...
    '$~N_>:' num2str(Nhigh) '$' newline ...
    'total counts : ' num2str(round(N1),'%.3g') ];
text(.98,.98,str,'units','normalized','fontsize',14,'horizontalalignment','right',...
    'verticalalignment','top','interpreter','latex');


Z2 = data.LatticeBin(2).Zbin;
Z2b = Z2;
Z2b(isnan(Z2)) = 0;
N2 = sum(Z2b,'all');

subplot(2,2,3);
imagesc(data.LatticeBin(2).n1,data.LatticeBin(2).n2,data.LatticeBin(2).Zbin);
colormap(cc);
axis equal tight
title('image 1');
caxis([0 2*Nthresh]);
title('image 2');
set(gca,'ydir','normal','box','on','linewidth',1);

ax_hB2=subplot(2,2,4);
pHistB12 = bar(1:100,1:100,'parent',ax_hB2,'linestyle','none',...
    'facecolor','k');
hold on
pHistBdivide2 = plot([1 1]*50,[0 100],'k-','parent',ax_hB2);
ylabel('occurences');
xlabel('counts per lattice site');
set(ax_hB2,'box','on','linewidth',.1,'fontsize',12,'units','normalized',...
    'XAxisLocation','bottom','YDir','normal','UserData','H1');
yyaxis right
pHistB22 = bar(1:100,1:100,'parent',ax_hB2,'linestyle','none',...
    'FaceColor',[0.6 0 0.5]);
ylabel('occurences');
ax_hB2.YColor=[0.6 0 0.5];   
x = data.LatticeHistogram(2).Centers;
xe = data.LatticeHistogram(2).Edges;
y = data.LatticeHistogram(2).N; 
xL = x<=Nthresh;
xH = ~xL;        
set(pHistB12,'XData',x(xL),'YData',y(xL));
set(pHistB22,'XData',x(xH),'YData',y(xH));        
pHistBdivide2.Parent.YAxis(1).Limits = [0 max(pHistB12.YData)*1.1];
pHistBdivide2.Parent.YAxis(2).Limits = [0 max(pHistB22.YData)*1.1];        
pHistBdivide2.Parent.XAxis.Limits = [0 max(x)*1.1];
set(pHistBdivide2,'Xdata',[1 1]*Nthresh,'Ydata',pHistBdivide2.Parent.YAxis(1).Limits);
title('image 2');
Nlow = sum(y(xL));
Nhigh = sum(y(xH));
str = ['threshold : ' num2str(Nthresh) newline ...
    '$N_>:' num2str(Nhigh) '$' newline ...
    'total counts : ' num2str(round(N2),'%.3g') ];
text(.98,.98,str,'units','normalized','fontsize',14,'horizontalalignment','right',...
    'verticalalignment','top','interpreter','latex');



end

