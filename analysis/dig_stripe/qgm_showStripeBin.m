function qgm_showStripeBin(qgmdata,xVar,opts)
    
Nseps = 7;
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
        % str = ['score=' num2str(scores(kk),'%.2e')];
        tScores(kk) = text(1,1,'','horizontalalignment','left',...
            'verticalalignment','middle','fontsize',10,...
            'color',myc);
    end
    
    % Text summary of stripe fit
    % str = ['$\lambda=' num2str(fout_s.L,'%.2f') ',' ...
    %     '\phi=2\pi\cdot' num2str(round(mod(fout_s.phi,2*pi)/(2*pi),2),'%.2f') ',' ...
    %     '\alpha = ' num2str(round(fout_s.B,2),'%.2f') '$'];
    tSummary = text(5,5,'','horizontalalignment','left',...
        'verticalalignment','bottom','fontsize',14,...
        'color',myc,'interpreter','latex','backgroundcolor',[0 0 0],'Margin',1,...
        'units','pixels');
    
    % Plot star for most in-focused plane
    pStar = plot(1,1,'pentagram','markersize',12,...
        'markerfacecolor',myc,'markeredgecolor',myc*.8);
    
    % %% Vertical Axis
    % ax2=subplot(5,5,[5 10 15 20]);
    % cla
    % 
    % % Data Plot
    % pZs = plot(Zs,ns,'k-','linewidth',1,'color','k');
    % hold on
    % nsFit = linspace(min(ns),max(ns),1e3);
    % 
    % % Fit Plot
    % pZsF = plot(feval(fout_s,nsFit),nsFit,'-','color',co(1,:),'linewidth',2);
    % set(gca,'fontsize',12,'YAxisLocation','right','YColor',co(1,:),'fontname','times',...
    %     'Xaxislocation','top')
    % ylabel('$n_2$ (site)','interpreter','latex');
    % ylim([min(ns) max(ns)])
    % plot(stripe_envelope(nsFit),nsFit,'--','color',co(1,:),'linewidth',1)
    % drawnow;
    % 
    % % Separations
    % for kk=1:Nseps
    %     pSeps2(kk) = plot(get(gca,'XLim'),[1 1]*seps(kk),'--','color',myc);
    % end
    % grid on
    % strR = ['$R^2 = ' num2str(gof_t.rsquare,'%.3f') '$'];
    % tRs = text(4,2,strR,'fontsize',10,'interpreter','latex','verticalalignment','bottom',...
    %     'units','pixels');
    % %% Horizontal Axis
    % ax3=subplot(5,5,[21 22 23 24]);
    % plot(nt,Zt,'k-','linewidth',1);
    % hold on
    % ntFit = linspace(min(nt),max(nt),1e3);
    % plot(ntFit,feval(fout_t,ntFit),'-','color',co(2,:),'linewidth',2);
    % set(gca,'ydir','normal','fontsize',12,'Xcolor',co(2,:),'fontname','times')
    % xlabel('$n_1$ (site)','interpreter','latex');
    % xlim([min(nt) max(nt)])
    % grid on
    % strR = ['$R^2 = ' num2str(gof_s.rsquare,'%.3f') '$'];
    % text(4,2,strR,'fontsize',10,'interpreter','latex','verticalalignment','bottom',...
    %     'units','pixels')
    % 
    % linkaxes([ax1 ax2],'y');
    % linkaxes([ax1 ax3],'x');
    % 
    % %% Score Axis
    % ax4 = subplot(5,5,25);
    % pScore = plot(centers,scores,'ko','markerfacecolor',[.5 .5 .5]);
    % xlabel('fringe position (site)');
    % ylabel('score');
    % set(gca,'fontsize',8);
    % yL =get(gca,'YLim');
    % ylim([0 yL(2)]);
    % xlim([min(n1) max(n1)])

    %% main Loop

    for nn = 1:length(qgmdata)
        n1 = [qgmdata(nn).LatticeBin.n1];
        n2 = [qgmdata(nn).LatticeBin.n2];
        Zb = [qgmdata(nn).LatticeBin.Zbin];
        BS = [qgmdata(nn).BinStripe];

        % Update Main Image
        set(hImg,'CData',Zb,'XData',n1,'YData',n2);

    
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
            str = ['score=' num2str(scores(kk),'%.2e')];
            set(tScores(kk),'String',str,'Position',[min(n1)+6 BS.Centers(kk)],'Visible','on')
        end
        set(pStar,'XData',min(n1)+1,'YData',BS.FocusCenter);


        % Update Summary 
        str = ['$\lambda=' num2str(BS.fitStripe.L,'%.2f') ',' ...
            '\phi=2\pi\cdot' num2str(round(mod(BS.fitStripe.phi,2*pi)/(2*pi),2),'%.2f') ',' ...
            '\alpha = ' num2str(round(BS.fitStripe.B,2),'%.2f') '$'];
        set(tSummary,'String',str);



    end
end

