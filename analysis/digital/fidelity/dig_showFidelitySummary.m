function [outputArg1,outputArg2] = dig_showFidelitySummary(fidelity)


digdata.Fidelity = fidelity;

%% Summary Figure

if length(digdata.FileNames)>1
    hF_out = figure;
    set(hF_out,'color','w','Name',['fidelity_summary']);
    clf
    hF_out.Position=[0 710 1300 250];
    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        t=uicontrol('style','text','string',opts.FigLabel,'fontsize',7,...
            'backgroundcolor','w','Position',[1 1 600 15],'horizontalalignment','left');
        t.Position(2) = t.Parent.Position(4)-t.Position(4)-2;
    end

    X= [digdata.X];
    co=get(gca,'colororder');
    subplot(141);
    p1=plot(X,[fidelity.N1],'ko','markerfacecolor',co(4,:),'color',...
        co(4,:)*.5,'linewidth',1);
    hold on
    p2=plot(X,[fidelity.N2],'^','markerfacecolor',co(5,:),'color',...
        co(5,:)*.5,'linewidth',1);
    ylabel('atom number')
    if isequal(digdata.xVar,'ExecutionDate')
        datetick x
    end
    legend([p1 p2],{'N1','N2'});
    yL=get(gca,'YLim');
    set(gca,'YLim',[0 yL(2)]);
    xlabel(digdata.xVar,'interpreter','none');

    subplot(142);
    pLost=plot(X,[fidelity.Nlost],'o','markerfacecolor',co(2,:),'color',...
        co(2,:)*.5,'linewidth',1);
    hold on
    pHop=plot(X,[fidelity.Nhop],'o','markerfacecolor',co(1,:),'color',...
        co(1,:)*.5,'linewidth',1);
    legend([pLost pHop],{'N lost','N hop'});
    ylabel('number')
    if isequal(digdata.xVar,'ExecutionDate')
        datetick x
    end
    yL=get(gca,'YLim');
    set(gca,'YLim',[0 yL(2)]);
    xlabel(digdata.xVar,'interpreter','none');

    subplot(143);
    pRLost=plot(X,[fidelity.Rlost],'o','markerfacecolor',co(2,:),'color',...
        co(2,:)*.5,'linewidth',1);
    hold on
    pRHop=plot(X,[fidelity.Rhop],'o','markerfacecolor',co(1,:),'color',...
        co(1,:)*.5,'linewidth',1);
    legend([pRLost pRHop],{'rate lost','rate hop'});
    ylabel('number')
    if isequal(digdata.xVar,'ExecutionDate')
        datetick x
    end
    yL=get(gca,'YLim');
    set(gca,'YLim',[0 yL(2)]);
    xlabel(digdata.xVar,'interpreter','none');

    subplot(144);
    dz = digdata.Zdig(:,:,:,1)-digdata.Zdig(:,:,:,2);
    defects = abs(dz);
    defects = sum(defects,3);
    imagesc(n1,n2,defects);
    colormap jet
    title('defect map');
    colorbar 
    axis equal tight

    xc = mean(digdata.Xc_site,'all');
    yc = mean(digdata.Yc_site,'all');

    xlim(xc+[-50 50])
    ylim(yc+[-50 50])


end

end

