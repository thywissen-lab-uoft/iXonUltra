function [hF] = pd_lattice_fluor(pd_data,opts)
%
hF=[];

hF = figure;
hF.Color = 'w';
hF.Name = 'photiode lattice fluor';
hF.Position = [100 100 750 500];
%% Analyze the lattice load traces
c=jet(length(pd_data));

xall = [];
yall = [];
zall = [];
for nn = 1:length(pd_data)      
    t = pd_data(nn).t;
    v = pd_data(nn).v;
    

    % Get the lattice photodiode information
    xlatt = v(:,7);
    ylatt = v(:,8);
    zlatt = v(:,9);
    
    i1 = find(xlatt>0.5*max(xlatt),1);    
    i2 = length(xlatt) - find(flip(xlatt)>0.5*max(xlatt),1) + 1;
    
    tlim = [t(i1)+10 t(i2)-10];
    
    ia = find(t>=(t(i1)+10),1);
    ib = find(t>=(t(i2)-10),1);

    xall = [xall; xlatt(ia:ib)];
    yall = [yall; ylatt(ia:ib)];
    zall = [zall; zlatt(ia:ib)];

    figure(hF);
    ax1=subplot(3,3,[1 2]);
    plot(t,xlatt,'-','linewidth',.5);
    hold on
    ylabel('xlatt - xlatt0 (mV)');   
    title('x lattice');
    xlim(tlim);
    
    ax2=subplot(3,3,[4 5]);
    plot(t,ylatt,'-','linewidth',.5);
    hold on
    ylabel('ylatt - ylatt0 (mV)');   
    title('y lattice');
    xlim(tlim);

    ax3=subplot(3,3,[7 8]);
    plot(t,zlatt,'-','linewidth',.5);
    hold on
    ylabel('zlatt - zlatt0 (mV)');   
    title('z lattice');
    xlim(tlim);
    xlabel('time (ms)');

end

subplot(3,3,3,'parent',hF);
pdx= fitdist(xall,'normal');
histfit(xall,50);
xlabel('x val (mV)');
set(gca,'YTick',[]);
str = [num2str(round(pdx.mu)) '\pm' num2str(round(pdx.sigma)) ' mV'];
text(.01,.98,str,'units','normalized','horizontalalignment','left',...
    'verticalalignment','top','fontsize',10);

subplot(3,3,6,'parent',hF);
pdy= fitdist(yall,'normal');
histfit(yall,50);
xlabel('y val (mV)');
set(gca,'YTick',[]);
str = [num2str(round(pdy.mu)) '\pm' num2str(round(pdy.sigma)) ' mV'];
text(.01,.98,str,'units','normalized','horizontalalignment','left',...
    'verticalalignment','top','fontsize',10);

subplot(3,3,9,'parent',hF);
pdz= fitdist(zall,'normal');
histfit(zall,50);
xlabel('z val (mV)');
set(gca,'YTick',[]);
str = [num2str(round(pdz.mu)) '\pm' num2str(round(pdz.sigma)) ' mV'];
text(.01,.98,str,'units','normalized','horizontalalignment','left',...
    'verticalalignment','top','fontsize',10);

if nargin >1 && isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)        
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left','parent',hF);
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=hF.Position(3);
    tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
end

end

