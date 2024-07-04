function [ixondata] = ixon_digitize(ixondata,threshold)

if ~isfield(ixondata,'LatticeBin')
    return;
end

% if size(ixondata(1).Z,3)

for kk=1:length(ixondata)
    LatticeDig = struct;
    for k = 1:length(ixondata(kk).LatticeBin)
 


        n1 = ixondata(kk).LatticeBin(k).n1;
        n2 = ixondata(kk).LatticeBin(k).n2;
        Zdig = ixondata(kk).LatticeBin(k).Zbin>=threshold; 
        Zb = ixondata(kk).LatticeBin(k).Zbin;
        Natoms = sum(sum(Zdig));        % Total number of atoms

        a1 =ixondata(kk).LatticeBin(k).a1;
        a2 =ixondata(kk).LatticeBin(k).a2;
        p1 =ixondata(kk).LatticeBin(k).p(1);
        p2 =ixondata(kk).LatticeBin(k).p(2);
        
        A = [a1 a2];                

        
        
        [nn1,nn2]=meshgrid(n1,n2);
        
        nn1=nn1(:);
        nn2=nn2(:);
        Zthere = Zdig(:);
        
        Counts=Zb(Zthere==1);
        Counts = Counts(:)';
        
        nn1(Zthere==0)=[];
        nn2(Zthere==0)=[];

        N=[nn1' ; nn2'];                      % All points
        p=[p1; p2];
        P = repmat(p,[1 length(N)]);        % Phase vector        
        Rn=A*(N+P);                                 % Positino of every atoms
        


        
      
          zY=sum(Zdig,2)';zY = zY/sum(zY);
        zX=sum(Zdig,1); zX = zX/sum(zX);

        % Calculate center of mass
        Xc_site=sum(zX.*n1);
        Yc_site=sum(zY.*n2);          

        % Calculate central second moment/variance and the standard
        % deviation
        X2_site=sum(zX.*(n1-Xc_site).^2); % x variance
        Xs_site=sqrt(X2_site); % standard deviation X
        Y2_site=sum(zY.*(n2-Yc_site).^2); % x variance
        Ys_site=sqrt(Y2_site); % standard deviation Y        
 


        [nn1,nn2]=meshgrid(n1,n2);
        X = (nn1+ixondata(kk).LatticeBin(k).p(1)).*a1(1) + ...
            (nn2+ixondata(kk).LatticeBin(k).p(2)).*a2(1);
        Y = (nn1+ixondata(kk).LatticeBin(k).p(1)).*a1(2) + ...
            (nn2+ixondata(kk).LatticeBin(k).p(2)).*a2(2);
        Xc_px = sum(X.*Zdig,'all')/sum(Zdig,'all');
        Yc_px = sum(Y.*Zdig,'all')/sum(Zdig,'all');
        X2_px = sum(X.^2.*Zdig,'all')/sum(Zdig,'all');
        Y2_px = sum(Y.^2.*Zdig,'all')/sum(Zdig,'all');
        Xs_px=sqrt(X2_px-Xc_px.^2); % standard deviation X
        Ys_px=sqrt(Y2_px-Yc_px.^2); % standard deviation Y        
        
        LatticeDig(k).Rn = Rn;
        LatticeDig(k).N = N;
        LatticeDig(k).Counts = Counts;
        LatticeDig(k).n1 = n1;
        LatticeDig(k).n2 = n2;
        LatticeDig(k).a1 = a1;
        LatticeDig(k).a2 = a2;
        LatticeDig(k).p1 = p1;
        LatticeDig(k).p2 = p2;

        LatticeDig(k).site2px = ixondata(kk).LatticeBin(k).site2px;
        LatticeDig(k).site2um = ixondata(kk).LatticeBin(k).site2um;
        LatticeDig(k).lattice_spacing_px =  ixondata(kk).LatticeBin(k).lattice_spacing_px;
        LatticeDig(k).lattice_spacing_um =  ixondata(kk).LatticeBin(k).lattice_spacing_um;
        LatticeDig(k).Zdig = Zdig;
        LatticeDig(k).Threshold = threshold;
        LatticeDig(k).Natoms = Natoms;

        LatticeDig(k).Xc_site = Xc_site;
        LatticeDig(k).Yc_site = Yc_site;
        LatticeDig(k).Xs_site = Xs_site;
        LatticeDig(k).Ys_site = Ys_site;   

        LatticeDig(k).Xc_px = Xc_px;
        LatticeDig(k).Yc_px = Yc_px;
        LatticeDig(k).Xs_px = Xs_px;
        LatticeDig(k).Ys_px = Ys_px;   

        LatticeDig(k).Xc_um = Xc_px*LatticeDig(k).lattice_spacing_um/LatticeDig(k).lattice_spacing_px;
        LatticeDig(k).Yc_um = Yc_px*LatticeDig(k).lattice_spacing_um/LatticeDig(k).lattice_spacing_px;
        LatticeDig(k).Xs_um = Xs_px*LatticeDig(k).lattice_spacing_um/LatticeDig(k).lattice_spacing_px;
        LatticeDig(k).Ys_um = Ys_px*LatticeDig(k).lattice_spacing_um/LatticeDig(k).lattice_spacing_px; 


    end
    ixondata(kk).LatticeDig = LatticeDig;
    
end

end

