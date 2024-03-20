function [pd_data,odt_summary,hF] = pd_odt_power(pd_data,opts)

hF=[];

hF = figure;
hF.Color = 'w';
hF.Name = 'photiode odt power';
hF.Position = [100 100 750 500];
%% Analyze the lattice load traces

odt1all = [];
odt2all = [];
for nn = 1:length(pd_data)      
    t = pd_data(nn).t;
    v = pd_data(nn).v;    

    % Get the odt photodiode information
    odt1 = v(:,3);
    odt2 = v(:,6);

    tlim = [-100 50];
    ia = find(t>=(tlim(1)),1);
    ib = find(t>=(tlim(2)),1);

    odt1all = [odt1all; odt1(ia:ib)];
    odt2all = [odt2all; odt2(ia:ib)];

    figure(hF);
    ax1=subplot(2,3,[1 2]);
    plot(t,odt1,'-','linewidth',.5);
    hold on
    ylabel('odt1 sucm (mV)');   
    title('odt1');
    xlim(tlim);
    
    ax2=subplot(2,3,[4 5]);
    plot(t,odt2,'-','linewidth',.5);
    hold on
    ylabel('odt2 sucm (mV)');   
    title('odt2');    
    xlim(tlim);
    xlabel('time (ms)');
    
    pd_data(nn).odt1_mV = mean(odt1(ia:ib));
    pd_data(nn).odt1_mV_err = std(odt1(ia:ib));
    pd_data(nn).odt2_mV = mean(odt2(ia:ib));
    pd_data(nn).odt2_mV_err = std(odt2(ia:ib));      
end
% 
subplot(2,3,3,'parent',hF);
pd1= fitdist(odt1all,'normal');
histfit(odt1all,50);
xlabel('odt1 val (mV)');
set(gca,'YTick',[]);
str = [num2str(round(pd1.mu,2)) '\pm' num2str(round(pd1.sigma,2)) ' mV'];
text(.01,.98,str,'units','normalized','horizontalalignment','left',...
    'verticalalignment','top','fontsize',10);

subplot(2,3,6,'parent',hF);
pd2= fitdist(odt2all,'normal');
histfit(odt2all,50);
xlabel('odt2 val (mV)');
set(gca,'YTick',[]);
str = [num2str(round(pd2.mu,2)) '\pm' num2str(round(pd2.sigma,2)) ' mV'];
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
%%

odt_summary = struct;
odt_summary.odt1_mV = pd1.mu;
odt_summary.odt1_mV_err = pd1.sigma;
odt_summary.odt2_mV = pd2.mu;
odt_summary.odt2_mV_err = pd2.sigma;

end

