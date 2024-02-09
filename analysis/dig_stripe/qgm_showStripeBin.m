function qgm_showStripeBin(qgmdata,xVar,opts)
    
Nseps = 10;
Ncenters = Nseps-1;

n1 = 1:50;
n2 = 1:50;
Zb = zeros(50,50);

%% Initialize Graphics
    hF=figure;
    co=get(gca,'colororder');
    myc = [255,140,0]/255;
    hF.Color='w';
    hF.Position(3:4) = [770 720];
    if (hF.Position(2)+hF.Position(4))>1000;hF.Position(2) = 100;end
    clf
    
    %% Image Axis
    ax1=subplot(5,5,[1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19]);

    % Image Plot
    hImg = imagesc(n1,n2,Zb);

    % Colorbar
    colormap([[0 0 0];winter; [1 0 0]]);    
    p = get(gca,'Position');
    c=colorbar('location','westoutside','fontsize',10,'color',[0 0 0 0],...
        'fontname','times');
    set(c.Label,'Color','k','String','counts/site')
    set(ax1,'position',p,'ydir','normal','fontsize',10,...
        'XAxisLocation','bottom','YColor',co(1,:),...
        'XColor',co(2,:),'fontname','times','yaxislocation','right')
    caxis(opts.ColorThreshold)
    hold on
    
    % Separation of Stripes
    for kk=1:Nseps
        pSeps(kk) = plot(1,1,'--','color',myc);
    end
    
    % Text label for each string
    for kk=1:Ncenters
        tScores(kk) = text(1,1,'','horizontalalignment','left',...
            'verticalalignment','middle','fontsize',10,...
            'color',myc);
    end
    
    % Text summary of stripe fit
    tSummary = text(5,5,'','horizontalalignment','left',...
        'verticalalignment','bottom','fontsize',14,...
        'color',myc,'interpreter','latex','backgroundcolor',[0 0 0],'Margin',1,...
        'units','pixels');
    
    % Plot star for most in-focused plane
    pStar = plot(1,1,'pentagram','markersize',12,...
        'markerfacecolor',myc,'markeredgecolor',myc*.8);
    
    %% Vertical Axis
    ax2=subplot(5,5,[5 10 15 20]);
    cla
    
    % Data Plot
    pZs = plot(1,1,'k-','linewidth',1,'color','k');
    hold on
    pZsF = plot(1,1,'-','color',co(1,:),'linewidth',2);
    set(gca,'fontsize',12,'YAxisLocation','right','YColor',co(1,:),'fontname','times',...
        'Xaxislocation','top')
    ylabel('$n_2$ (site)','interpreter','latex');
    ylim([0 1])
    
%     pSEnv = plot(1,1,'--','color',co(1,:),'linewidth',1);
%     drawnow;
%     
    % Separations
    for kk=1:Nseps
        pSeps2(kk) = plot(1,1,'--','color',myc);
    end
    grid on
    tRs = text(4,2,'','fontsize',10,'interpreter','latex','verticalalignment','bottom',...
        'units','pixels');
    
    
    % %% Horizontal Axis
    ax3=subplot(5,5,[21 22 23 24]);
    pZt = plot(1,1,'k-','linewidth',1);
    hold on
    pZtF = plot(1,1,'-','color',co(2,:),'linewidth',2);
    set(gca,'ydir','normal','fontsize',12,'Xcolor',co(2,:),'fontname','times')
    xlabel('$n_1$ (site)','interpreter','latex');
    xlim([0 1])
    grid on
    tRt = text(4,2,'','fontsize',10,'interpreter','latex','verticalalignment','bottom',...
        'units','pixels');

    % %% Score Axis
    ax4 = subplot(5,5,25);
    pScore = plot(1,1,'ko','markerfacecolor',[.5 .5 .5]);
    xlabel('fringe position (site)');
    ylabel('score');
    set(gca,'fontsize',8);
    ylim([0 1]);
    xlim([0 1])

    %% main Loop

    for nn = 1:length(qgmdata)
        n1 = [qgmdata(nn).LatticeBin.n1];
        n2 = [qgmdata(nn).LatticeBin.n2];
        Zb = [qgmdata(nn).LatticeBin.Zbin];
        Zb2 = Zb;
        Zb2(Zb<opts.ColorThreshold(1)) = 0;
        Zb2(isnan(Zb)) = 0;
        
        BS = [qgmdata(nn).BinStripe];

        % Update Main Image
        set(hImg,'CData',Zb,'XData',n1,'YData',n2);
        set(ax1,'XLim',[min(n1) max(n1)],'YLim',[min(n2) max(n2)]);
    
        % Update Separation Bars
        for kk = 1:length(pSeps)
            pSeps(kk).Visible='off';
        end
        for kk=1:length([BS.Separations])
            set(pSeps(kk),'XData',[min(n1) max(n1)],'YData',[1 1]*BS.Separations(kk),'Visible','on');
        end

        % Update Score Text
        for kk = 1:length(tScores)
            tScores(kk).Visible='off';
        end
        for kk=1:length([BS.Centers])
            str = ['score=' num2str(BS.Scores(kk),'%.2e')];
            set(tScores(kk),'String',str,'Position',[min(n1)+6 BS.Centers(kk)],'Visible','on')
        end
        set(pStar,'XData',min(n1)+1,'YData',BS.FocusCenter);

        % Update Summary 
        str = ['$\lambda=' num2str(BS.FitStripe.L,'%.2f') ',' ...
            '\phi=2\pi\cdot' num2str(mod(BS.FitStripe.phi,2*pi)/(2*pi),'%.2f') ',' ...
            '\alpha = ' num2str(BS.FitStripe.B,'%.2f') '$'];
        set(tSummary,'String',str);
        
        set(pZs,'XData',sum(Zb2,2)','YData',n2);
        set(pZsF,'XData',feval(BS.FitStripe,n2)','YData',n2);
        set(ax2,'YLim',[min(n2) max(n2)]);
        
        % Update Separation Bars
        for kk = 1:length(pSeps2)
            pSeps2(kk).Visible='off';
        end
        for kk=1:length([BS.Separations])
            set(pSeps2(kk),'XData',get(ax2,'XLim'),'YData',[1 1]*BS.Separations(kk),'Visible','on');
        end
        set(tRs,'String',['$R^2 = ' num2str(BS.RSquareStripe,'%.3f') '$']);

        set(pZt,'YData',sum(Zb2,1)','XData',n1);
        set(pZtF,'YData',feval(BS.FitTransverse,n1)','XData',n1);
        set(ax3,'XLim',[min(n1) max(n1)]);
        set(tRt,'String',['$R^2 = ' num2str(BS.RSquareTransverse,'%.3f') '$']);


        set(pScore,'XData',[BS.Centers],'YData',[BS.Scores]);
        try
        set(ax4,'XLim',[min(n1) max(n1)],'YLim',[0 max([BS.Scores])*1.2]);
        end
        drawnow;
        
        
    end
end

