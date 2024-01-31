function data = ixon_binnedHistogram(data,Nbins)
        if ~isfield(data,'LatticeDig')
           return;
        end      
       
        if nargin == 1
           Nbins = 200;
        end
%         Nbins = histBtbl.Data(1,2);    
        
        Zall = data.LatticeBin(1).Zbin;
        Zall = Zall(:);
        Zall(Zall==0) = NaN;                     
        [N,edges] = histcounts(Zall,Nbins);        
        centers = (edges(1:end-1) + edges(2:end))/2;
        LatticeHistogram = struct;
        LatticeHistogram.Edges = edges;
        LatticeHistogram.Centers = centers;
        LatticeHistogram.N = N;  
        data.LatticeHistogram(1) = LatticeHistogram; 
        for kk=2:length(data.LatticeBin)
            Zall = data.LatticeBin(kk).Zbin;
            Zall = Zall(:);
%           Zall(Zall==0)=[];            
            [N,edges] = histcounts(Zall,data.LatticeHistogram(1).Edges);
            centers = (edges(1:end-1) + edges(2:end))/2;
            LatticeHistogram = struct;
            LatticeHistogram.Edges = edges;
            LatticeHistogram.Centers = centers;
            LatticeHistogram.N = N;  
            data.LatticeHistogram(kk) = LatticeHistogram;            
        end  
        
end

