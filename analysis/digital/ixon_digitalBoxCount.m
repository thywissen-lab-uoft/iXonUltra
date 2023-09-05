function ixondata = ixon_digitalBoxCount(ixondata)
    fprintf('Performing digitized box count analysis ...');    

%     data.LatticeDig(kk) = ixon_digitalBoxCount(data.LatticeDig(kk));

    for kk=1:length(ixondata)
        LatticeBin = ixondata(kk).LatticeDig;
        LatticeDig = struct;
        for k=1:length(LatticeBin)
            Zdig = LatticeBin(k).Zdig;      
            
            x = LatticeBin(k).n1;
            y = LatticeBin(k).n2;
            
            Natoms = sum(sum(Zdig));        % Total number of atoms
            Nsites = length(Zdig(:));       % Total number of sites        
%             filling = Natoms/Nsites;        % Average filling fraction
            
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

            LatticeDig(k).Natoms = Natoms;
%             LatticeDig(k).Filling = filling;
            LatticeDig(k).Xc = Xc;
            LatticeDig(k).Yc = Yc;
            LatticeDig(k).Xs = Xs;
            LatticeDig(k).Ys = Ys;   
        end         
        ixondata(kk).LatticeDig = LatticeDig;
    end
    disp('done');
end

