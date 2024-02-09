function [hF] = qgm_binnedTotalHistogram(qgmdata,opts)
    if ~isfield(qgmdata,'LatticeBin')
       return;
    end   
    
%% Options
    if nargin == 1
        opts = struct;
    end
    
    if ~isfield(opts,'ROI')
        opts.ROI = 'max';
    end   
   
    if ~isfield(opts,'Nthresh')
       opts.Nthresh = 1000; 
    end
    
    if ~isfield(opts,'Bins')
        opts.Bins = 100;
    end   
    
    if ~isfield(opts,'doAnimate')
       opts.doAnimate = 0; 
    end
%% Get ROI
    
    n1 = qgmdata(1).LatticeBin(1).n1;
    n2 = qgmdata(1).LatticeBin(1).n2;
    
    if isequal(opts.ROI,'max')
       opts.ROI = [min(n1) max(n1) min(n2) max(n2)]; 
    end    
    R = opts.ROI;
    %% Prepare Data
    Zall = zeros(length(n2),length(n1));
    for nn = 1:length(qgmdata)        
        Zthis = qgmdata(nn).LatticeBin(1).Zbin;    
        Zall(:,:,nn) =  Zthis;
    end
    
    in1i = find(n1==R(1),1);in1f = find(n1==R(2),1);    
    in2i = find(n2==R(3),1);in2f = find(n2==R(4),1);
    
     
  [N,edges] = histcounts(Zall(in2i:in2f,in1i:in1f,:),opts.Bins);  
    centers = (edges(1:end-1) + edges(2:end))/2;           
    iL = centers<=opts.Nthresh;
    iH = ~iL; 
    
    %% Initialize Graphics
    hF = figure;
    hF.Color='w';
    hF.Position= [100 100 1300 400];
    hF.Name = 'Binned Histogram';
    
    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        tFig=uicontrol('style','text','string',opts.FigLabel,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left');
        tFig.Position(4)=tFig.Extent(4);
        tFig.Position(3)=hF.Position(3);
        tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
    end    
    
    % Histogram Axis
    ax1 = subplot(121);        
    
    % Low Counts
    pHistB1 = bar(centers(iL),N(iL),'linestyle','none',...
        'facecolor','k');
    xlim([200 max(edges)]);    
    ylabel('occurences');
    xlabel('counts per lattice site');
    hold on
    
    % High Counts
    yyaxis right
    pHistB2 = bar(centers(iH),N(iH),'linestyle','none',...
    'FaceColor',[0.6 0 0.5]);
    set(gca,'box','on','YColor',[0.6 0 0.5]*.8);
    ylabel('occurences');
    
    % String Labels
    str1 = ['ROI : [' num2str(R(1)) ' ' num2str(R(2)) ' ' ...
        num2str(R(3)) ' ' num2str(R(4)) ']' ...
        newline num2str(size(Zall,3)) ' images' ];
    t1 = text(.98,.98,str1,'units','normalized','verticalalignment',...
        'top','horizontalalignment','right','fontsize',8);

    % Image
    ax2 = subplot(122);
    Zc = Zall;
    Zc(isnan(Zc)) = 0;
    hImg = imagesc(n1,n2,mean(Zc,3));    
    rectangle('Position',[R(1) R(3) R(2)-R(1) R(4)-R(3)],...
        'EdgeColor','r')
    xlabel('n1 sites');
    ylabel('n2 sites');
    c=colorbar;
    c.Label.String = 'counts/site';
    axis equal tight
    set(gca,'box','on','linewidth',1,'fontsize',10,'ydir','normal');
    ca = [0 0 0];       
    cb = [0.7 .1 .6];
    cc = [linspace(ca(1),cb(1),1000)' ...
        linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
    colormap(hF,cc);
    
    %% Populate with data
    
    if opts.doAnimate
        for nn = 1:length(qgmdata)
            [N,~] = histcounts(Zall(in2i:in2f,in1i:in1f,nn),opts.Bins);  
            set(pHistB1,'Xdata',centers(iL),'Ydata',N(iL));
            set(pHistB2,'Xdata',centers(iH),'Ydata',N(iH));
            set(hImg,'Cdata',Zall(:,:,nn));
               str1 = ['ROI : [' num2str(R(1)) ' ' num2str(R(2)) ' ' ...
        num2str(R(3)) ' ' num2str(R(4)) ']' ...
        newline 'image ' num2str(nn) '/' num2str(length(qgmdata))];
    
            set(t1,'String',str1);
            
            frame=getframe(hF);
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);              
            filename = fullfile(opts.saveDir,opts.filename);            
            switch nn
                case 1
                    imwrite(A,map,filename,'gif','LoopCount',...
                        Inf,'DelayTime',1);
                case length(qgmdata)
                    imwrite(A,map,filename,'gif','WriteMode',...
                        'append','DelayTime',1);
                otherwise
                    imwrite(A,map,filename,'gif','WriteMode',...
                        'append','DelayTime',.1);
             end
               
        
            
        end

        
    end
    
    %% Last Image
    
                [N,~] = histcounts(Zall(in2i:in2f,in1i:in1f,:),opts.Bins);  
            set(pHistB1,'Xdata',centers(iL),'Ydata',N(iL));
            set(pHistB2,'Xdata',centers(iH),'Ydata',N(iH));
            Zc = Zall;
            Zc(isnan(Zc)) = 0;
            c.Label.String = 'average counts/site';

            set(hImg,'Cdata',mean(Zc,3));
    str1 = ['ROI : [' num2str(R(1)) ' ' num2str(R(2)) ' ' ...
        num2str(R(3)) ' ' num2str(R(4)) ']' ...
        newline num2str(size(Zall,3)) ' images' ];

    
            set(t1,'String',str1);
    
end

