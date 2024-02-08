% ixon_QGM.m
%
% Author : CF Fujiwara
%
% This script is the primary analysis code for ixon images which are single
% plane.

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));   

qgmdata = struct;

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
end

%% Find the Lattice
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

%% Plot Lattice Analysis

if ixon_doQGM_FindLattice
    [LatticeK,hF] = ixon_showLatticeK(ixondata);   
end

%% 

reassignBadK = 1;
useAverageK = 0;

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

%% 
if ixon_doQGM_FindLattice
    
    if ~useAverageK
        [hF] = ixon_showLatticeA(ixondata);
    end
    
    
    
end

%% Assign Lattice Data to QGM data

%% Bin Data
if ixon_doQGM_Bin
    for n=1:length(ixondata)
        for kk=1:size(ixondata(n).Z,3)
            opts = struct;
            a1 = ixondata(n).LatticePhase(kk).a1;
            a2 = ixondata(n).LatticePhase(kk).a2;                        
            p1 = ixondata(n).LatticePhase(kk).p1;
            p2 = ixondata(n).LatticePhase(kk).p2;        
            opts.ScaleFactor = 3;    
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
        end    
        qgmdata(n).LatticeBin = ixondata(n).LatticeBin;
    end
end

%% Bin Stripe
if ixon_doQGM_BinStripe
    LGuess = 25;
    ColorThreshold = [1000 3000];
    
    if ~isfield(ixondata,'LatticeBin')
        return;
    end
    clear out
    for n = 1:length(ixondata)                
        for kk = 1:length(ixondata(n).LatticeBin)
            n1 = ixondata(n).LatticeBin(kk).n1;
            n2 = ixondata(n).LatticeBin(kk).n2;
            Zb = ixondata(n).LatticeBin(kk).Zbin;    
            opts_stripe.LGuess = LGuess;
            opts_stripe.FigNum=3000+10*(n-1)+kk-1;
                        opts_stripe.FigNum=3000;

            opts_stripe.ColorThreshold = ColorThreshold;
            
            [out(kk),hF_bin_stripe] = ixon_fitStripe_dig(n1,n2,Zb,opts_stripe);
%                   exportgraphics(gcf,'testAnimated.gif','Append',true);
            frame=getframe(hF_bin_stripe);
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);  
            
            filename = fullfile(saveDir,'binstripe_75.gif');
            
            if out(kk).ModDepth>=.75 && out(kk).Counts>0.5e6
                switch n
                    case 1
                        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
                    case length(ixondata)
                        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
                    otherwise
                        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.1);
                end
            end

        end
        ixondata(n).BinStripe = out;    
        qgmdata(n).BinStripe = out;
    end
end  

%% Digitization Stuff
if ixon_doQGM_Digitize
    ixon_main_digital;
end