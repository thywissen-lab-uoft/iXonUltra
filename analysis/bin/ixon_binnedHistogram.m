function [data,Zbb] = ixon_binnedHistogram(data,Nbins)
        if ~isfield(data,'LatticeBin')
           return;
        end      
       
        if nargin == 1
           Nbins = 200;
        end
%         Nbins = histBtbl.Data(1,2);    
% Zbb = {};
for nn = 1:length(data)        
        Zall = data(nn).LatticeBin(1).Zbin;
        Zall = Zall(:);
%         Zbb{end+1} = Zall;

        Zall(Zall==0) = NaN;                     
        [N,edges] = histcounts(Zall,Nbins);        
        centers = (edges(1:end-1) + edges(2:end))/2;
        LatticeHistogram = struct;
        LatticeHistogram.Edges = edges;
        LatticeHistogram.Centers = centers;
        LatticeHistogram.N = N;  
        data(nn).LatticeHistogram(1) = LatticeHistogram; 
        for kk=2:length(data(nn).LatticeBin)
            Zall = data(nn).LatticeBin(kk).Zbin;
            Zall = Zall(:);
%             Zbb{end+1} = Zall;

%             Zbb = [Zbb; Zall];
%           Zall(Zall==0)=[];            
            [N,edges] = histcounts(Zall,data(nn).LatticeHistogram(1).Edges);
            centers = (edges(1:end-1) + edges(2:end))/2;
            LatticeHistogram = struct;
            LatticeHistogram.Edges = edges;
            LatticeHistogram.Centers = centers;
            LatticeHistogram.N = N;  
            data(nn).LatticeHistogram(kk) = LatticeHistogram;            
        end  
end

end

