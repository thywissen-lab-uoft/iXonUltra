function [digdata] = bin_makeDigDataScaled(bindata,opts)
    P = [bindata.Params];
    F = [bindata.Flags];
    U = [bindata.Units];
    B = [bindata.LatticeBin];
    PH =[bindata.LatticePhase];
    
    Zdig = zeros(size(bindata(1).LatticeBin(1).Zbin,1),...
        size(bindata(1).LatticeBin(1).Zbin,2),length(bindata));
    n1 = bindata(1).LatticeBin.n1;
    n2 = bindata(1).LatticeBin.n2;
    
    thresh = [];
    
    thresh_med     = median([bindata.ScaledThreshold]);
    thresh_abs_med = median([bindata.ScaledThreshold].*[bindata.ScaledCentroid]);
    
    if opts.DigAve
            thresh     = thresh_med;
            thresh_abs = thresh_abs_med;
            
        Threshtype = 'CompensatedAve';
    else
        Threshtype = 'CompensatedInd';
    end


    for nn=1:length(bindata)
        Zscaled = bindata(nn).LatticeBin(1).ZbinScaled;
        
        if opts.DigAve
            Zdig(:,:,nn) = Zscaled>=thresh;  
        else
            thresh(nn)=bindata(nn).ScaledThreshold;
            thresh_abs(nn) = thresh(nn)*bindata(nn).ScaledCentroid;
            if thresh_abs(nn)<2000
                thresh(nn) = thresh_med;
                thresh_abs(nn) = thresh_abs_med;
            end
            Zdig(:,:,nn) = Zscaled>=thresh(nn);  

        end
        
        Natoms(nn) = sum(Zdig(:,:,nn),[1 2]);
    end
    
    
    digdata                     = struct;    
    digdata.SourceDirectory     = unique({bindata.SourceDirectory});
    digdata.FileNames           = {bindata.Name}';
    digdata.ScaledThreshold     = thresh;
    digdata.Threshold           = thresh_abs;
    digdata.ThresholdingType    = Threshtype;
    digdata.xVar                = opts.xVar;
    digdata.X                   = [P.(opts.xVar)];
    digdata.Params              = P;
    digdata.Units               = U;
    digdata.Flags               = F;      
    
    % Actual Digital Data
    digdata.Zdig                = logical(Zdig);
    digdata.n1                  = n1;
    digdata.n2                  = n2;
    digdata.Natoms              = Natoms;       % Number of atoms 
        
    % Basis Information
    digdata.a1                  = [B.a1];       % Lattice Vector 1
    digdata.a2                  = [B.a2];       % Lattice Vector 2    
    digdata.p1                  = [PH.p1];      % Lattice Phase 1
    digdata.p2                  = [PH.p2];      % Lattice Phase 2
    
    [digdata] = dig_basic(digdata);
end

