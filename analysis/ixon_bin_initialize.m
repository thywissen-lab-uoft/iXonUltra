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
    BC = [ixondata.BoxCount];
    N  = [BC.Ncounts];

    for nn=1:length(ixondata)
        fprintf(['(' num2str(nn) '/' num2str(numel(ixondata))...
            ') lattice phase']);        
        k1 = ixondata(nn).LatticeK(1).k1;
        k2 = ixondata(nn).LatticeK(1).k2;     

        if N(nn)<0.25*max(N)
            LatticeK.BadLattice(nn) = 1;
        end

        if (ixon_doQGM_reassignBadK && LatticeK.BadLattice(nn)) || ixon_doQGM_useAverageK
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
    bin_options = struct;    
    ixondata = ixon_ProcessImagesToBin(ixondata);    
end  
%% Initialize bindata    
bindata = ixon_makeBinData(ixondata,saveOpts);
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
    try if ~exist(saveOpts.saveDir,'dir');mkdir(saveOpts.saveDir);end;end
    filename = fullfile(saveOpts.saveDir,'bindata.mat');
    disp(['Saving ' filename ' ...']);
    save(filename, 'bindata');
end