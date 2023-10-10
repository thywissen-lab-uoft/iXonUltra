qgmdata = struct;

%% Options
qgm_doFindLattice = 1;
qgm_doBinLattice = 1;
qgm_doDigitize   = 1;
qgm_DigitizationThreshhold = 1500;
qgm_doDigitalAnalysis = 1;

%%
if qgm_doFindLattice
    for n=1:length(ixondata)
        for kk=1:size(ixondata(n).Zf,3)
            tic;
            fprintf(['(' num2str(kk) '/' num2str(size(ixondata(n).Zf,3)) ') Fitting reciprocal lattice ...']);
            ixondata(n).LatticeK(kk) = findLatticeK(ixondata(n).f,ixondata(n).f,ixondata(n).Zf(:,:,kk),ixondata(n).ProcessOptions);                
            k1 = ixondata(n).LatticeK(kk).k1;
            k2 = ixondata(n).LatticeK(kk).k2;          
            t2 = toc;
            fprintf([' done (' num2str(t2,2) ' sec.)' ' phase ...']);
            tic;
            ixondata(n).LatticePhase(kk) = findLatticePhase(ixondata(n).X,ixondata(n).Y,ixondata(n).Z,k1,k2);              
            t=toc;
            disp([' done (' num2str(t,2) ' sec.)']);       
      
            qgmdata(n).ProcessOptions = ixondata(n).ProcessOptions;
            qgmdata(n).Z = ixondata(n).Z;
            qgmdata(n).LatticeK(kk) = ixondata(n).LatticeK(kk);
            qgmdata(n).LatticePhase(kk) = ixondata(n).LatticePhase(kk);
            qgmdata(n).Params = ixondata(n).Params;
            qgmdata(n).Units = ixondata(n).Units;
            qgmdata(n).Flags = ixondata(n).Flags;
        end
    end    
end

%% 
if qgm_doBinLattice
    for n=1:length(ixondata)
        for kk=1:size(ixondata(n).Z,3)
            opts = struct;
            a1 = ixondata(n).LatticePhase(kk).a1;
            a2 = ixondata(n).LatticePhase(kk).a2;                        
            p1 = ixondata(n).LatticePhase(kk).p1;
            p2 = ixondata(n).LatticePhase(kk).p2;        
            opts.ScaleFactor = 4;    
            opts.a1 = a1;
            opts.a2 = a2;
            opts.p1 = p1;
            opts.p2 = p2;     
            if isfield(ixondata(n),'RotationMask')
               opts.Mask =  ixondata(n).RotationMask;
            end
            ROI=ixondata(n).ROI;
            ix_1 = find(ixondata(n).X>=ROI(1),1);
            ix_2 = find(ixondata(n).X>=ROI(2),1);
            iy_1 = find(ixondata(n).Y>=ROI(3),1);
            iy_2 = find(ixondata(n).Y>=ROI(4),1);
            x = ixondata(n).X(ix_1:ix_2);
            y = ixondata(n).Y(iy_1:iy_2);   
            z = ixondata(n).Z(iy_1:iy_2,ix_1:ix_2,kk);    
            tic;
            fprintf(['(' num2str(kk) '/' num2str(size(ixondata(n).Zf,3)) ...
                ') binning into lattice ...']);   
            ixondata(n).LatticeBin(kk) = binLattice(x,y,z,opts); 
            t2=toc;
            disp(['done (' num2str(t2,3) ' sec.)']);                
            qgmdata(n).LatticeBin(kk) = ixondata(n).LatticeBin(kk);
        end    
    end
end

%%

if qgm_doDigitize
   [ixondata] = ixon_digitize(ixondata,qgm_DigitizationThreshhold);
   for n=1:length(ixondata)
       qgmdata(n).LatticeDig = ixondata(n).LatticeDig;
   end
end
            
%% 
if qgm_doDigitalAnalysis
   qgm_DigitalAnalysis; 
end