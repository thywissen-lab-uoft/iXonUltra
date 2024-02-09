function [hF] = qgm_binnedTotalHistogram(qgmdata,Nbins)
    if ~isfield(qgmdata,'LatticeBin')
       return;
    end      
    Nthresh=1000;
    if nargin == 1
       Nbins = 200;
    end

    Zall = [];
    for nn = 1:length(qgmdata)        
        Zthis = qgmdata(nn).LatticeBin(1).Zbin;    
        
%         Zthis = Zthis(80:120,100:150);
        Zall = [Zall; Zthis(:)];
    end
    
    [N,edges] = histcounts(Zall,Nbins);  
    centers = (edges(1:end-1) + edges(2:end))/2;           
    iL = centers<=Nthresh;
    iH = ~iL;          
                
    hF = figure;
    hF.Color='w';
        
    pHistB1 = bar(centers(iL),N(iL),'linestyle','none',...
        'facecolor','k');
    xlim([200 max(edges)]);    
    ylabel('occurences');
    xlabel('counts per lattice site');
    hold on
    
    plot([1 1]*Nthresh,get(gca,'YLim'),'k-');
    yyaxis right
    pHistB2 = bar(centers(iH),N(iH),'linestyle','none',...
    'FaceColor',[0.6 0 0.5]);
    set(gca,'box','on','YColor',[0.6 0 0.5]*.8);
    ylabel('occurences');

end

