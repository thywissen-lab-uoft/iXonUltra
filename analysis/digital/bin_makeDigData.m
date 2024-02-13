function [digdata] = bin_makeDigData(bindata,opts)
    P = [bindata.Params];
    F = [bindata.Flags];
    U = [bindata.Units];
    D = [bindata.LatticeDig];

    
    Zdig = zeros(size(bindata(1).LatticeBin.Zbin,1),...
        size(bindata(1).LatticeBin.Zbin,2),length(bindata));
    n1 = bindata(1).LatticeBin.n1;
    n2 = bindata(1).LatticeBin.n2;

    for nn=1:length(bindata)
        Zdig(:,:,nn) = [bindata(nn).LatticeDig.Zdig];        
    end
    
    LD = [bindata.LatticeDig];
    
    digdata = struct;    
    digdata.SourceDirectory = unique({bindata.SourceDirectory});
    digdata.FileNames = {bindata.Name}';
    digdata.Threshold = [LD.Threshold];
    digdata.xVar = opts.xVar;
    digdata.X = [P.(opts.xVar)];
    digdata.Params = P;
    digdata.Units = U;
    digdata.Flags = F;      
    digdata.Zdig = logical(Zdig);
    digdata.n1 = n1;
    digdata.n2 = n2;
    digdata.Natoms = [D.Natoms];
    digdata.Xc = [D.Xc];
    digdata.Yc = [D.Yc];
    digdata.Xs = [D.Xs];
    digdata.Ys = [D.Ys];
end

