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

    digdata.lattice_spacing_px = [D.lattice_spacing_px];
    digdata.lattice_spacing_um = [D.lattice_spacing_um];
    digdata.a1 = [D.a1];
    digdata.a2 = [D.a2];
    digdata.p1 = [D.p1];
    digdata.p2 = [D.p2];

    digdata.site2px = {D.a1};
    digdata.site2um = {D.a2};

    digdata.Xc_site = [D.Xc_site];
    digdata.Yc_site = [D.Yc_site];
    digdata.Xs_site = [D.Xs_site];
    digdata.Ys_site = [D.Ys_site];

    digdata.Xc_px = [D.Xc_px];
    digdata.Yc_px = [D.Yc_px];
    digdata.Xs_px = [D.Xs_px];
    digdata.Ys_px = [D.Ys_px];

    digdata.Xc_um = [D.Xc_um];
    digdata.Yc_um = [D.Yc_um];
    digdata.Xs_um = [D.Xs_um];
    digdata.Ys_um = [D.Ys_um];
end

