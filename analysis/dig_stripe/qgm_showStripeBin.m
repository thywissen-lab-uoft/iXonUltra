function qgm_showStripeBin(qgmdata,xVar,opts)
    

%% Initialize Graphics
    hF=figure;
    co=get(gca,'colororder');
    myc = [255,140,0]/255;
    hF.Color='w';
    hF.Position(3:4) = [770 720];
    if (hF.Position(2)+hF.Position(4))>1000;hF.Position(2) = 100;end
    clf
    
    % Image Plot
    ax1=subplot(5,5,[1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19]);
    hImg = imagesc(n1,n2,Zb);
    % axis equal tight
    colormap([[0 0 0];winter; [1 0 0]]);
    
    p = get(gca,'Position');
    c=colorbar('location','westoutside','fontsize',10,'color',[0 0 0 0],...
        'fontname','times');
    c.Label.Color='k';
    c.Label.String ='counts/site';
    set(gca,'position',p)
    set(gca,'ydir','normal','fontsize',10,'XAxisLocation','bottom','YColor',co(1,:),...
        'XColor',co(2,:),'fontname','times','yaxislocation','right')
    caxis(opts.ColorThreshold)
    hold on
    
    % Lines for indicating the absence of atoms
    for kk=1:length(seps)
        pSeps(kk) = plot(get(gca,'XLim'),[1 1]*seps(kk),'--','color',myc);
    end
    
    % Text label for each string
    for kk=1:length(centers)
        str = ['score=' num2str(scores(kk),'%.2e')];
        tScores(kk) = text(min(nt)+6,centers(kk),str,'horizontalalignment','left',...
            'verticalalignment','middle','fontsize',10,...
            'color',myc);
    end
    
    % Text summary of stripe fit
    str = ['$\lambda=' num2str(fout_s.L,'%.2f') ',' ...
        '\phi=2\pi\cdot' num2str(round(mod(fout_s.phi,2*pi)/(2*pi),2),'%.2f') ',' ...
        '\alpha = ' num2str(round(fout_s.B,2),'%.2f') '$'];
    tSummary = text(5,5,str,'horizontalalignment','left',...
        'verticalalignment','bottom','fontsize',14,...
        'color',myc,'interpreter','latex','backgroundcolor',[0 0 0],'Margin',1,...
        'units','pixels')
    
    % Plot star for most in-focused plane
    pStar = plot(min(nt)+2,focus_center,'pentagram','markersize',12,...
        'markerfacecolor',myc,'markeredgecolor',myc*.8);
    
    ax2=subplot(5,5,[5 10 15 20]);
    cla
    pZs = plot(Zs,ns,'k-','linewidth',1,'color','k');
    hold on
    nsFit = linspace(min(ns),max(ns),1e3);
    pZsF = plot(feval(fout_s,nsFit),nsFit,'-','color',co(1,:),'linewidth',2);
    set(gca,'fontsize',12,'YAxisLocation','right','YColor',co(1,:),'fontname','times',...
        'Xaxislocation','top')
    ylabel('$n_2$ (site)','interpreter','latex');
    ylim([min(ns) max(ns)])
    plot(stripe_envelope(nsFit),nsFit,'--','color',co(1,:),'linewidth',1)
    drawnow;
    for kk=1:length(seps)
        plot(get(gca,'XLim'),[1 1]*seps(kk),'--','color',myc)
    end
    grid on
    strR = ['$R^2 = ' num2str(gof_t.rsquare,'%.3f') '$'];
    text(4,2,strR,'fontsize',10,'interpreter','latex','verticalalignment','bottom',...
        'units','pixels');
    
    ax3=subplot(5,5,[21 22 23 24]);
    plot(nt,Zt,'k-','linewidth',1);
    hold on
    ntFit = linspace(min(nt),max(nt),1e3);
    plot(ntFit,feval(fout_t,ntFit),'-','color',co(2,:),'linewidth',2);
    set(gca,'ydir','normal','fontsize',12,'Xcolor',co(2,:),'fontname','times')
    xlabel('$n_1$ (site)','interpreter','latex');
    xlim([min(nt) max(nt)])
    grid on
    strR = ['$R^2 = ' num2str(gof_s.rsquare,'%.3f') '$'];
    text(4,2,strR,'fontsize',10,'interpreter','latex','verticalalignment','bottom',...
        'units','pixels')
    
    linkaxes([ax1 ax2],'y');
    linkaxes([ax1 ax3],'x');
    
    subplot(5,5,25);
    plot(centers,scores,'ko','markerfacecolor',[.5 .5 .5]);
    xlabel('fringe position (site)');
    ylabel('score');
    set(gca,'fontsize',8);
    yL =get(gca,'YLim');
    ylim([0 yL(2)]);
    xlim([min(n1) max(n1)])
end

