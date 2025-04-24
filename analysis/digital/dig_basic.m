function [digdata] = dig_basic(digdata)
    if ~isfield(digdata,'Lattice_px') || ~isfield(digdata,'Lattice_um')
        digdata.Lattice_px          = 2.68;
        digdata.Lattice_um          = 1.054/2;
    end
    a_px = digdata.Lattice_px;
    a_um = digdata.Lattice_um;  
%     keyboard
%     %Some legacy issue
%     if length(digdata.Ratom(:,1)) == 1
%         digdata.Ratom = digdata.Ratom';
%     end
%     if length(digdata.Natoms(:,1)) == 1
%         digdata.Natoms = digdata.Natoms';
%     end
    
    
    for nn=1:length(digdata.FileNames)
        for rr=1:size(digdata.Zdig,4)
            Ratom = digdata.Ratom{nn,rr};
            
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
              
            digdata.Xc_px(nn,rr) = Xc_px;
            digdata.Xc_site(nn,rr) = Xc_site;
            digdata.Xc_um(nn,rr) = Xc_um;
            digdata.Xs_px(nn,rr) = Xs_px;
            digdata.Xs_site(nn,rr) = Xs_site;
            digdata.Xs_um(nn,rr) = Xs_um;
            
            digdata.Yc_px(nn,rr) = Yc_px;
            digdata.Yc_site(nn,rr) = Yc_site;
            digdata.Yc_um(nn,rr) = Yc_um;
            digdata.Ys_px(nn,rr) = Ys_px;
            digdata.Ys_site(nn,rr) = Ys_site;
            digdata.Ys_um(nn,rr) = Ys_um;
    
            digdata.nPeakGauss(nn,rr) = digdata.Natoms(nn,rr)./(sqrt(2*pi*Xs_site.^2).*sqrt(2*pi*Ys_site.^2));
        end
    end  
end

