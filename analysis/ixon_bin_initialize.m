% ixon_bin_initialize.m
%
% Author : CF Fujiwara
%
% This script is the primary analysis code for ixon images which are single
% plane.

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));  

% Initialize outputs
qgmdata = struct;

% Some flags
ixon_doQGM_FindLattice = 1;
reassignBadK = 1;
useAverageK = 0;
ixon_doQGM_Bin = 1;

%% Initial checks
% Check to make sure that PSF deconvolution had been done 
if ~exist('ixondata')
   error(['This code assume you have run ixon_main.m first, which should ' ...
       'define the variable ixondata.']);
end
img_opts = [ixondata.ProcessOptions];
for kk=1:length(img_opts)
    if ~img_opts(kk).doPSF
        error(['Run ' num2str(kk) ' of ' length(img_opts) ' was not ' ...
            'sharpened with the PSF deconvlution.']);
    end
end

%% Find Lattice K
if ixon_doQGM_FindLattice
    for n = 1:length(ixondata)
        tic
        fprintf(['(' num2str(n) '/' num2str(numel(ixondata))...
                ') lattice wavevector']);
        for kk=1:size(ixondata(n).Zf,3)
            fprintf(['...' num2str(kk)]);
            ixondata(n).LatticeK(kk) = findLatticeK(...
                ixondata(n).f,ixondata(n).f,...
                ixondata(n).Zf(:,:,kk)); 
        end 
        disp([' done (' num2str(toc,'%.2f') 's)']);
    end
end

[LatticeK,hF_LatticeK] = ixon_showLatticeK(ixondata);   


%% Find Lattice Phase

if ixon_doQGM_FindLattice
    for nn=1:length(ixondata)
        fprintf(['(' num2str(nn) '/' num2str(numel(ixondata))...
            ') lattice phase']);        
        k1 = ixondata(nn).LatticeK(1).k1;
        k2 = ixondata(nn).LatticeK(1).k2;          
        if (reassignBadK && LatticeK.BadLattice(nn)) || useAverageK
            k1 = LatticeK.k1;
            k2 = LatticeK.k2;
        end        
        % Make sure column vector
        k1 = reshape(k1,[2 1]);
        k2 = reshape(k2,[2 1]);
        tic
        for kk=1:size(ixondata(nn).Zf,3)
            fprintf(['...' num2str(kk)]);
            ixondata(nn).LatticePhase(kk) = findLatticePhase(...
                ixondata(nn).X,ixondata(nn).Y,ixondata(nn).Z,k1,k2);              
        end
        disp([' done (' num2str(toc,'%.2f') 's)']);        
    end
end

%% Bin Data
if ixon_doQGM_Bin
    if isfield(ixondata,'LatticeBin')
    ixondata = rmfield(ixondata,'LatticeBin');
    end
    for n=1:length(ixondata)
        fprintf(['(' num2str(n) '/' num2str(numel(ixondata))...
            ') lattice binning']);
        tic
        for kk=1:size(ixondata(n).Z,3)
            fprintf(['...' num2str(kk)]);

            opts = struct;
            a1 = ixondata(n).LatticePhase(kk).a1;
            a2 = ixondata(n).LatticePhase(kk).a2;                        
            p1 = ixondata(n).LatticePhase(kk).p1;
            p2 = ixondata(n).LatticePhase(kk).p2;        
            opts.ScaleFactor = 2;    
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
            ixondata(n).LatticeBin(kk) = binLattice(x,y,z,opts); 
        end    
        disp([' done (' num2str(toc,'%.2f') 's)']);        
    end
end  

%% Initialize QGM Data    
for nn = 1:length(ixondata)
    qgmdata(nn).Date = ixondata(nn).Date;
    qgmdata(nn).Name = ixondata(nn).Name;
    qgmdata(nn).Params = ixondata(nn).Params;
    qgmdata(nn).Units = ixondata(nn).Units;
    qgmdata(nn).Flags = ixondata(nn).Flags;
    qgmdata(nn).CameraInformation = ixondata(nn).CameraInformation;
    qgmdata(nn).AcquisitionInformation = ixondata(nn).AcquisitionInformation;
    qgmdata(nn).AcquisitionDescription = ixondata(nn).AcquisitionDescription;
    qgmdata(nn).ROI = ixondata(nn).ROI;
    qgmdata(nn).ProcessOptions = ixondata(nn).ProcessOptions;  
    qgmdata(nn).LatticeK = ixondata(nn).LatticeK;
    qgmdata(nn).LatticePhase = ixondata(nn).LatticePhase;
    qgmdata(nn).LatticeBin = ixondata(nn).LatticeBin;
end

%% Save Figures

hF_LatticeVectors = ixon_showLatticeA(ixondata);
hF_LatticePhase = ixon_showLatticePhase(ixondata);    
[LatticeK,hF_LatticeK] = ixon_showLatticeK(ixondata);   

if ixon_doSave        
    ixon_saveFigure2(hF_LatticeK,...
        'ixon_LatticeK',saveOpts);
    ixon_saveFigure2(hF_LatticeVectors,...
        'ixon_LatticeVectors',saveOpts);
    ixon_saveFigure2(hF_LatticePhase,...
        'ixon_LatticePhase',saveOpts);      
end

%% Save QGM Data
if ixon_doSave       
    filename = fullfile(saveOpts.saveDir,'qgmdata.mat');
    save(filename, 'qgmdata');
end