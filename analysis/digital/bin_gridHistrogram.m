function hF = bin_gridHistrogram(bindata,opts)
    if ~isfield(bindata,'LatticeBin')
       return;
    end   
    
      if ~isfield(opts,'Nthresh')
       opts.Nthresh = 2000; 
    end
    
    if ~isfield(opts,'Bins')
        opts.Bins = 100;
    end 
end

