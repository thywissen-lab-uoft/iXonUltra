function hF = bin_gridHistogram(bindata,opts)
    if ~isfield(bindata,'LatticeBin')
       return;
    end   
    
    if nargin ==1
        opts=struct;
    end
    
    if ~isfield(opts,'Nthresh')
       opts.Nthresh = 2000; 
    end
    
    if ~isfield(opts,'NumBins')
        opts.NumBins = 80;
    end 
    
   if ~isfield(opts,'NumGrid')
      opts.NumGrid=[3 3]; 
   end
    
    edges = linspace(0,6*opts.Nthresh,opts.NumBins);
    n1 = bindata(1).LatticeBin(1).n1;
    n2 = bindata(1).LatticeBin(1).n2;
    %% Get All Data
    
    Zall = zeros(length(n2),length(n1),length(bindata));
    for nn = 1:length(bindata)        
        Zthis = bindata(nn).LatticeBin(1).Zbin;    
        Zall(:,:,nn) =  Zthis;
    end
    
    Zallcopy = Zall;
    Zallcopy(isnan(Zallcopy))=0;
    
    Zsum = sum(Zallcopy,3);
    
    ZsumX = sum(Zsum,1);
    
    
    
    n1c = sum(n1.*ZsumX)/sum(ZsumX);
    x2 = sum(n1.^2.*ZsumX)/sum(ZsumX);    
    sx = sqrt(x2-n1c^2);    
    n1i = floor(n1c-2*sx)-2;
    n1f = ceil(n1c+2*sx)+2;
    n1_inspect = n1i:1:n1f;

    ZsumY = sum(Zsum,2)';
    n2c = sum(n2.*ZsumY)/sum(ZsumY);
    y2 = sum(n2.^2.*ZsumY)/sum(ZsumY);    
    sy = sqrt(y2-n2c^2);     
    n2i = floor(n2c-2*sy)-2;
    n2f = ceil(n2c+2*sy)+2;    
    n2_inspect = n2i:1:n2f;  
    
    
    %% Prepare Data
   
    dNc = floor(length(n1_inspect)/opts.NumGrid(1));
    dNr = floor(length(n2_inspect)/opts.NumGrid(2));

    
    N={};
%     edges={};
    R={};
    for rr=1:opts.NumGrid(1)
        for cc = 1:opts.NumGrid(2)
            
            n10 = n1_inspect(1)+(cc-1)*dNc;
            n1f = n1_inspect(1)+cc*dNc;
            
            if cc==opts.NumGrid(1)
                n1f = n1_inspect(end);
            end
            
            n20 = n2_inspect(1)+(rr-1)*dNr;
            n2f = n2_inspect(1)+rr*dNr;
            
            if rr==opts.NumGrid(2)
                n2f = n2_inspect(end);
            end
            
            
       
            
            in1i = find(n1==n10,1);in1f = find(n1==n1f,1);    
            in2i = find(n2==n20,1);in2f = find(n2==n2f,1);
            
            R{end+1}=[n10 n20 n1f-n10 n2f-n20];
            
            
            [N{rr,cc}] = histcounts(Zall(in2i:in2f,in1i:in1f,:),edges);           
          
       
        end
    end
    
    
         centers = (edges(1:end-1) + edges(2:end))/2;           
            iL = centers<=opts.Nthresh;
            iH = ~iL; 
        %% Initialize Graphics
    hF = figure;
    hF.Color='w';
    hF.Position= [50 50 1400 400];
    hF.Name = 'Binned Histogram';
    

    p1 = uipanel('parent',hF,'units','normalized','position',[0 0 .7 1],...
        'backgroundcolor','w');
    
    % Histogram Axis
    j=1;
    chsv=hsv(numel(N));
    
    for rr=1:size(N,1)
        for cc=1:size(N,2)
            ind = (size(N,1)-rr)*size(N,1) + cc;
            subplot(opts.NumGrid(1),opts.NumGrid(1),ind,'parent',p1);    
            
            Nme = N{rr,cc};
            Nme(1)=NaN;
            % Low Counts
            pHistB1 = bar(centers(iL),Nme(iL),'linestyle','none',...
                'facecolor','k');
            xlim([5 max(edges)]);    
            ylabel('occurences');
            xlabel('counts per lattice site');
            hold on

            % High Counts
            yyaxis right
            pHistB2 = bar(centers(iH),Nme(iH),'linestyle','none',...
            'FaceColor',[0.6 0 0.5]);
            set(gca,'box','on','YColor',[0.6 0 0.5]*.8);
            ylabel('occurences');
                title(['(' num2str(rr) ',' num2str(cc) ')'],'color',chsv(j,:));

            j=j+1;
        end
    end
    
        if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        tFig=uicontrol('style','text','string',opts.FigLabel,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left','parent',p1);
        tFig.Position(4)=tFig.Extent(4);
        tFig.Position(3)=hF.Position(3);
        tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
    end    
    
    
    %% Image
    p2 = uipanel('parent',hF,'units','normalized','position',[.7 0 .3 1],...
        'backgroundcolor','w');
    axes('parent',p2);
    imagesc(n1,n2,mean(Zall,3));

    
 xlabel('n1 sites');
    ylabel('n2 sites');
    c=colorbar;
    c.Label.String = 'counts/site';
    caxis([0 2*opts.Nthresh]);
    
    axis equal tight
    set(gca,'box','on','linewidth',1,'fontsize',10,'ydir','normal');
    ca = [0 0 0];       
    cb = [0.7 .1 .6];
    cc = [linspace(ca(1),cb(1),1000)' ...
        linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
    colormap(hF,cc);
    hold on
    
    for nn=1:length(R)
       rectangle('Position',R{nn},'EdgeColor',chsv(nn,:)); 
    end
    title('average binned imaged');
    

end

