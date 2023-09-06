function [ixondata] = ixon_digitize(ixondata,threshold)

if ~isfield(ixondata,'LatticeBin')
    return;
end
for kk=1:length(ixondata)
    LatticeDig = struct;
    for k = 1:length(ixondata.LatticeBin)
        x = ixondata(kk).LatticeBin(k).n1;
        y = ixondata(kk).LatticeBin(k).n2;
        Zdig = ixondata(kk).LatticeBin(k).Zbin>=threshold; 
        Natoms = sum(sum(Zdig));        % Total number of atoms




        zY=sum(Zdig,2)';zY = zY/sum(zY);
        zX=sum(Zdig,1); zX = zX/sum(zX);

        % Calculate center of mass
        Xc=sum(zX.*x);
        Yc=sum(zY.*y);          

        % Calculate central second moment/variance and the standard
        % deviation
        X2=sum(zX.*(x-Xc).^2); % x variance
        Xs=sqrt(X2); % standard deviation X
        Y2=sum(zY.*(y-Yc).^2); % x variance
        Ys=sqrt(Y2); % standard deviation Y               

        LatticeDig.Zdig = Zdig;
        LatticeDig.Natoms = Natoms;
        LatticeDig.n1 = x;
        LatticeDig.n2 = y;
        LatticeDig.Xc = Xc;
        LatticeDig.Yc = Yc;
        LatticeDig.Xs = Xs;
        LatticeDig.Ys = Ys;   
    end
    ixondata(kk).LatticeDig = LatticeDig;
    
end

end

