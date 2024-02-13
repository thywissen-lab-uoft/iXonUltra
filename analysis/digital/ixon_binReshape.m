function bindata = ixon_binReshape(bindata)

% INCOMPLETE
% Get all bounds on the lattice sites indeces
n1i = min(LatticeDig(1).n1);n1f = max(LatticeDig(1).n1);
n2i = min(LatticeDig(1).n2);n2f = max(LatticeDig(1).n2);

for n = 1:length(bindata)
    
end
% Get all bounds on the lattice sites indeces
n1i = min(LatticeDig(1).n1);n1f = max(LatticeDig(1).n1);
n2i = min(LatticeDig(1).n2);n2f = max(LatticeDig(1).n2);

% Find the maximu
for n =1:length(LatticeDig)
    n1i = min([LatticeDig(n).n1 n1i]);
    n1f = max([LatticeDig(n).n1 n1f]);
    n2i = min([LatticeDig(n).n2 n2i]);
    n2f = max([LatticeDig(n).n2 n2f]);
end

% Redfine the lattice vectors
n1 = n1i:n1f;
n2 = n2i:n2f;

% Iterate through all images and center the digitzal in the new lattice
% sites
Zdig_all = zeros(length(n2),length(n1),length(LatticeDig));
for n = 1:length(LatticeDig)
    n1i0 = LatticeDig(n).n1(1); 
    n2i0 = LatticeDig(n).n2(1); 
    dI1 = length(LatticeDig(n).n1);
    dI2 = length(LatticeDig(n).n2);    
    i1 = find(n1==n1i0,1);
    i2 = find(n2==n2i0,1);
    Zthis = zeros(length(n2),length(n1));
    Zthis(i2:(i2+dI2-1),i1:(i1+dI1-1)) = LatticeDig(n).Zdig;
    
    Zdig_all(:,:,n) = Zthis;
end

digdata.n1 = n1;
digdata.n2 = n2;
digdata.Z = Zdig_all;

end

