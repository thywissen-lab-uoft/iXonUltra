function [bindata] = bin_recenter(bindata)

LatticeBin = [bindata.LatticeBin];


% Get all bounds on the lattice sites indeces
n1i = min(LatticeBin(1).n1);n1f = max(LatticeBin(1).n1);
n2i = min(LatticeBin(1).n2);n2f = max(LatticeBin(1).n2);

% Find the maximum
for n =1:length(LatticeBin)
    n1i = min([LatticeBin(n).n1 n1i]);
    n1f = max([LatticeBin(n).n1 n1f]);
    n2i = min([LatticeBin(n).n2 n2i]);
    n2f = max([LatticeBin(n).n2 n2f]);
end

% Redfine the lattice vectors
n1 = n1i:n1f;
n2 = n2i:n2f;

% Iterate through all images and center the digitzal in the new lattice
% sites
Zbin_all = zeros(length(n2),length(n1),length(bindata),length(bindata(1).LatticeBin));
for n = 1:length(bindata)
    LB = bindata(n).LatticeBin;    
    for kk=1:length(LB)
        n1i0 = LB(kk).n1(1); 
        n2i0 = LB(kk).n2(1); 
        dI1 = length(LB(kk).n1);
        dI2 = length(LB(kk).n2);    
        i1 = find(n1==n1i0,1);
        i2 = find(n2==n2i0,1);
        Zthis = nan(length(n2),length(n1));
        Zthis(i2:(i2+dI2-1),i1:(i1+dI1-1)) = LB(kk).ZbinRaw;    
        Zbin_all(:,:,n,kk) = Zthis;
    end
end

for n = 1:length(bindata)
    K = length(bindata(n).LatticeBin);
    for kk=1:K
    bindata(n).LatticeBin(kk).n1 = n1; 
    bindata(n).LatticeBin(kk).n2 = n2; 
    bindata(n).LatticeBin(kk).ZbinRaw = Zbin_all(:,:,n,kk);
    bindata(n).LatticeBin(kk).Zbin = Zbin_all(:,:,n,kk);
    end
end


end

