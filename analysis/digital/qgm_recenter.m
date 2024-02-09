function [qgmdata] = qgm_recenter(qgmdata)

LatticeBin = [qgmdata.LatticeBin];

% Get all bounds on the lattice sites indeces
n1i = min(LatticeBin(1).n1);n1f = max(LatticeBin(1).n1);
n2i = min(LatticeBin(1).n2);n2f = max(LatticeBin(1).n2);

% Find the maximu
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
Zbin_all = zeros(length(n2),length(n1),length(LatticeBin));
for n = 1:length(LatticeBin)
    n1i0 = LatticeBin(n).n1(1); 
    n2i0 = LatticeBin(n).n2(1); 
    dI1 = length(LatticeBin(n).n1);
    dI2 = length(LatticeBin(n).n2);    
    i1 = find(n1==n1i0,1);
    i2 = find(n2==n2i0,1);
    Zthis = nan(length(n2),length(n1));
    
    Zthis(i2:(i2+dI2-1),i1:(i1+dI1-1)) = LatticeBin(n).Zbin;
    
    Zbin_all(:,:,n) = Zthis;
end

for n = 1:length(qgmdata)
    qgmdata(n).LatticeBin.n1 = n1; 
    qgmdata(n).LatticeBin.n2 = n2; 
    qgmdata(n).LatticeBin.Zbin = Zbin_all(:,:,n);
end


end

