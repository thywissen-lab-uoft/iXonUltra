function [digdata] = dig_basic(digdata)
    if ~isfield(digdata,'Lattice_px') || ~isfield(digdata,'Lattice_um')
        digdata.Lattice_px          = 2.68;
        digdata.Lattice_um          = 1.054/2;
    end
    a_px = digdata.Lattice_px;
    a_um = digdata.Lattice_um;
    


            
    for nn=1:length(digdata.FileNames)
        Ratom = digdata.Ratom{nn};
        
        Xc_px   = mean(Ratom(1,:));
        Xc_site = Xc_px/a_px;
        Xc_um   = Xc_site*a_um;
        
        Xs_px   = std(Ratom(1,:));
        Xs_site = Xs_px/a_px;
        Xs_um   = Xs_site*a_um;
        
        Yc_px   = mean(Ratom(2,:));
        Yc_site = Yc_px/a_px;
        Yc_um   = Yc_site*a_um;
        
        Ys_px   = std(Ratom(2,:));
        Ys_site = Ys_px/a_px;
        Ys_um   = Ys_site*a_um;        
          
        digdata.Xc_px(nn) = Xc_px;
        digdata.Xc_site(nn) = Xc_site;
        digdata.Xc_um(nn) = Xc_um;
        digdata.Xs_px(nn) = Xs_px;
        digdata.Xs_site(nn) = Xs_site;
        digdata.Xs_um(nn) = Xs_um;
        
        digdata.Yc_px(nn) = Yc_px;
        digdata.Yc_site(nn) = Yc_site;
        digdata.Yc_um(nn) = Yc_um;
        digdata.Ys_px(nn) = Ys_px;
        digdata.Ys_site(nn) = Ys_site;
        digdata.Ys_um(nn) = Ys_um;

        digdata.nPeakGauss(nn) = digdata.Natoms(nn)./(sqrt(2*pi*Xs_site.^2).*sqrt(2*pi*Ys_site.^2));
    end  
end

