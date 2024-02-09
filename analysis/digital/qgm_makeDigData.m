function [digdata] = qgm_makeDigData(qgmdata,opts)
    P = [qgmdata.Params];
    F = [qgmdata.Flags];
    U = [qgmdata.Units];
    D = [qgmdata.LatticeDig];

    
    Zdig = zeros(size(qgmdata(1).LatticeBin.Zbin,1),...
        size(qgmdata(1).LatticeBin.Zbin,2),length(qgmdata));
    n1 = qgmdata(1).LatticeBin.n1;
    n2 = qgmdata(1).LatticeBin.n2;

    for nn=1:length(qgmdata)
        Zdig(:,:,nn) = [qgmdata(nn).LatticeDig.Zdig];        
    end
    
    digdata = struct;
    digdata.FileNames = {qgmdata.Name}';
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

