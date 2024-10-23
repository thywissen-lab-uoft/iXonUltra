function bin_showStripeBin(bindata,xVar,opts)
    
Nseps = 10;
Ncenters = Nseps-1;

n1 = 1:50;
n2 = 1:50;
Zb = zeros(50,50);

if nargin <3
   opts=struct;
end

if nargin<2 || isempty(xVar)
    xVar = 'ExecutionDate';
end

if ~isfield(opts,'Threshold')
    opts.Threshold = [1000 5000];
end

if ~isfield(opts,'FigLabel')
    opts.FigLabel = [];
end

if ~isfield(opts,'FigureNumber')
    opts.FigureNumber = 30001;
end

if ~isfield(opts,'doSave')
   opts.doSave = 0; 
end

if ~isfield(opts,'SumIndex')
    opts.SumIndex = 2;
end


%% Initialize Graphics
    hF=figure(opts.FigureNumber);
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
    caxis(opts.Threshold)
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
    
    % Text summary of stripe fit
    tVar = text(.99,.99,'','horizontalalignment','right',...
        'verticalalignment','top','fontsize',14,...
        'color',myc,'interpreter','latex','backgroundcolor',[0 0 0],'Margin',1,...
        'units','normalized');  
    
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
    
    if opts.SumIndex == 2
        ylabel('$n_2$ (site)','interpreter','latex');
    else
        ylabel('$n_1$ (site)','interpreter','latex');
    end
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
    
    if opts.SumIndex == 2
    xlabel('$n_1$ (site)','interpreter','latex');
    else
    xlabel('$n_2$ (site)','interpreter','latex');
    end
        
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
    
    
    if ~isempty(opts.FigLabel)
    
        t=uicontrol('style','text','string',opts.FigLabel,'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left');
        t.Position(4)=t.Extent(4);
        t.Position(3)=hF.Position(3);
        t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
    end
    %% Main Loop

    for nn = 1:length(bindata)
        
        if opts.SumIndex == 2
        n1 = [bindata(nn).LatticeBin(1).n1];
        n2 = [bindata(nn).LatticeBin(1).n2];
        Zb = [bindata(nn).LatticeBin(1).Zbin];
        else
        n2 = [bindata(nn).LatticeBin(1).n1];
        n1 = [bindata(nn).LatticeBin(1).n2];
        Zb = [bindata(nn).LatticeBin(1).Zbin]';
        end
        
        
        Zb2 = Zb;
        Zb2(Zb<opts.Threshold(1)) = 0;
        Zb2(isnan(Zb)) = 0;
        
        BS = [bindata(nn).BinStripe(1)];

        % Update Main Image
        set(hImg,'CData',Zb,'XData',n1,'YData',n2);
        set(ax1,'XLim',[min(n1) max(n1)],'YLim',[min(n2) max(n2)]);
    
        % Update Separation Bars
        for kk = 1:length(pSeps)
            pSeps(kk).Visible='off';
        end
        for kk=1:length([BS(1).Separations])
            set(pSeps(kk),'XData',[min(n1) max(n1)],'YData',[1 1]*BS(1).Separations(kk),'Visible','on');
        end

        % Update Score Text
        for kk = 1:length(tScores)
            tScores(kk).Visible='off';
        end
        for kk=1:length([BS(1).Centers])
            str = ['score=' num2str(BS(1).Scores(kk),'%.2e')];
            set(tScores(kk),'String',str,'Position',[min(n1)+6 BS(1).Centers(kk)],'Visible','on')
        end
        set(pStar,'XData',min(n1)+1,'YData',BS(1).FocusCenter);

        % Update Summary 
        str = ['$\lambda=' num2str(BS(1).FitStripe.L,'%.2f') ',' ...
            '\phi=2\pi\cdot' num2str(mod(BS(1).FitStripe.phi,2*pi)/(2*pi),'%.2f') ',' ...
            '\alpha = ' num2str(BS(1).FitStripe.B,'%.2f') '$'];
        set(tSummary,'String',str);
        
        % Label image with date/xVar in top right (for analysis only)s
        if isequal(xVar,'ExecutionDate')
            strVar=[xVar ': ' datestr(bindata(nn).Params.(xVar),'YYYY-mm-DD_HH-MM-SS')];          % Variable string
        else
            strVar = [xVar ':' num2str(bindata(nn).Params.(xVar)) ];
        end
        set(tVar,'String',strVar,'interpreter','none');
        set(pZs,'XData',sum(Zb2,2)','YData',n2);
        set(pZsF,'XData',feval(BS(1).FitStripe,n2)','YData',n2);
        set(ax2,'YLim',[min(n2) max(n2)]);
        
        % Update Separation Bars
        for kk = 1:length(pSeps2)
            pSeps2(kk).Visible='off';
        end
        for kk=1:length([BS(1).Separations])
            set(pSeps2(kk),'XData',get(ax2,'XLim'),'YData',[1 1]*BS(1).Separations(kk),'Visible','on');
        end
        
        
        set(tRs,'String',['$R^2 = ' num2str(BS(1).RSquareStripe,'%.3f') '$']);

        set(pZt,'YData',sum(Zb2,1)','XData',n1);
        set(pZtF,'YData',feval(BS(1).FitTransverse,n1)','XData',n1);
        set(ax3,'XLim',[min(n1) max(n1)]);
        set(tRt,'String',['$R^2 = ' num2str(BS(1).RSquareTransverse,'%.3f') '$']);


        set(pScore,'XData',[BS(1).Centers],'YData',[BS(1).Scores]);
        try
        set(ax4,'XLim',[min(n1) max(n1)],'YLim',[0 max([BS(1).Scores])*1.2]);
        end
        drawnow;
        
        
        if opts.doSave
            frame=getframe(hF);
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);              
            filename = fullfile(opts.saveDir,opts.filename);            
            switch nn
                case 1
                    imwrite(A,map,filename,'gif','LoopCount',...
                        Inf,'DelayTime',1);
                case length(bindata)
                    imwrite(A,map,filename,'gif','WriteMode',...
                        'append','DelayTime',1);
                otherwise
                    imwrite(A,map,filename,'gif','WriteMode',...
                        'append','DelayTime',.1);
            end
        end
        
        
    end
end

