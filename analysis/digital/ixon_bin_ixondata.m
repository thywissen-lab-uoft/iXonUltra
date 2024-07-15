function ixondata = ixon_bin_ixondata(ixondata,opts)

if nargin==1
   opts=struct; 
end

if ~isfield(opts,'LatticeSource')
   opts.LatticeBasis = 'auto'; 
end

if  ~isfield(opts,'ROI')
   ROI = [1 512 1 512]; 
end

binOpts = struct;
binOpts.ROI = ROI;

for nn=1:length(ixondata)
    for kk=1:size(ixondata(nn).Z,3)
        X = ixondata(nn).X;
        Y = ixondata(nn).Y;
        Z = ixondata(nn).Z;
        
        if isequal(opts.LatticeSource,'auto')
            binOpts.a1 = ixondata(nn).LatticePhase(kk).a1; 
            binOpts.a2 = ixondata(nn).LatticePhase(kk).a2; 
            binOpts.p1 = ixondata(nn).LatticePhase(kk).p1; 
            binOpts.p2 = ixondata(nn).LatticePhase(kk).p2; 
        end
          
        ix_1 = find(X>=ROI(1),1);
        ix_2 = find(X>=ROI(2),1);
        iy_1 = find(Y>=ROI(3),1);
        iy_2 = find(Y>=ROI(4),1);
        
        x = X(ix_1:ix_2);
        y = Y(iy_1:iy_2);   
        z = Z(iy_1:iy_2,ix_1:ix_2,kk);   
        tic;
        fprintf(['(' num2str(kk) '/' num2str(size(data.Zf,3)) ...
            ') binning into lattice ...']);    

        data.LatticeBin(kk) = binLattice(x,y,z,binOpts); 
        t2=toc;
        disp(['done (' num2str(t2,3) ' sec.)']);
    end     
end
        
        
end

