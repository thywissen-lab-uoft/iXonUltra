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
            'sharpeend with the PSF deconvlution.']);
    end
end

%% Initialize the structure
qgmdata = struct;
%% Options

% These are the analysis options, typically these all need to be on in
% order to digitze the image
qgm_doFindLattice = 1;          
qgm_doBinLattice = 1;



qgm_doDigitalAnalysis = 1;

%% Using the FFT find the k vectors of the lattice and the phase
% Using the FFT, find the kvectors of the lattice and fit them.  The
% k-vectors give the lattice spacing.
%
% Using the k-vectors caluate the phase of the phase by fully calcuating
% the Fourier Transform. (FFT algorithm isn't precise enough)
%
% If the SNR is poor, the k-vectors fit may fail. This is determined by
% comparing the measured k-vectors to a predetermined "close enough" one
% which is set by the magnification, lattice angle, and lattice beams.
% These dont change unless the optics are altered.
%
% If the fitted k-vector differ from these "close enough" values too much
% the data pont is assigned the average fitted kvector at the end.
if qgm_doFindLattice
    
% This flag checks whether the fitted k-vectors are close to an
% approximation
doApproxCheck = 1;
badLatticeAction = 'Use Average';

% Define some expected k-vectors "close enough"
k10 = [0.3749; 0];k20 = [0.0031; 0.3743];    

% Average fitted k-vectors
k1Avg = [0;0];
k2Avg = [0;0];

% Keep track of which indices have a good k-vector fit
goodLatticeInds = ones(1,length(ixondata));
    
disp(repmat('-',1,60));disp(repmat('-',1,60)); 
disp(' Finding the lattice vectors and phase ...');
for n=1:length(ixondata)
    fprintf([num2str(n,'%02.f') '/' num2str(length(ixondata),'%02.f') ' ']);
    for kk=1:size(ixondata(n).Zf,3)
        tic;
        fprintf(['' num2str(kk) '/' num2str(size(ixondata(n).Zf,3)) ...
            ' ']);
        % Fit lattice wavevectors from FFT
        ixondata(n).LatticeK(kk) = findLatticeK(ixondata(n).f,...
            ixondata(n).f,...
            ixondata(n).Zf(:,:,kk),...
            ixondata(n).ProcessOptions);   
        
        % Display the result
        k1 = ixondata(n).LatticeK(kk).k1;k2 = ixondata(n).LatticeK(kk).k2;
        fprintf(['k1=(' num2str(round(k1(1),3)) ...
            ',' num2str(round(k1(2),3)) ') ']);
        fprintf(['k2=(' num2str(round(k2(1),3)) ...
            ',' num2str(round(k2(2),3)) ') | ']);
        
        errNorm1 = abs(norm(k1)-norm(k10))/norm(k10);
        errNorm2 = abs(norm(k2)-norm(k20))/norm(k20);
        
        dtheta1 = acos((sum(k1.*k10))/(norm(k1)*norm(k10)))*180/pi;
        dtheta2 = acos((sum(k2.*k20))/(norm(k2)*norm(k20)))*180/pi;

        if errNorm1>0.1 
           warning(['Fitted k1 has norm deviation of ' num2str(errNorm1*100,'%02.1f') '%.']);
           goodLatticeInds(n) = 0;
        end
        if errNorm2>0.1 
           warning(['Ftted k1 has norm deviation of ' num2str(errNorm2*100,'%02.1f') '%.']);
         goodLatticeInds(n) = 0;
        end
        if dtheta1>5
            warning(['k1 has angle deviation of ' num2str(dtheta1,'%02.1f') ' deg']);
             goodLatticeInds(n) = 0;
        end
        if dtheta2>5
            warning(['k2 has angle deviation of ' num2str(dtheta2,'%02.1f') ' deg']);
            goodLatticeInds(n) = 0;
        end            
        
        if goodLatticeInds(n)
           k1Avg = k1Avg + k1; 
           k2Avg = k2Avg + k2; 
        end

        % Calculate lattice phase from Fourier Transform
        ixondata(n).LatticePhase(kk) = findLatticePhase(...
            ixondata(n).X,ixondata(n).Y,ixondata(n).Z,...
            k1,k2);   
        
        % Display the result
        a1 = ixondata(n).LatticePhase(kk).a1;a2 = ixondata(n).LatticePhase(kk).a2;
        p1 = ixondata(n).LatticePhase(kk).p1;p2 = ixondata(n).LatticePhase(kk).p2;        
        fprintf(['a1=(' num2str(round(a1(1),2)) ...
            ',' num2str(round(a1(2),2)) ')']);
        fprintf(['a2=(' num2str(round(a2(1),2)) ...
            ',' num2str(round(a1(1),2)) ')']);
        fprintf(['phase=2*pi*(' num2str(round(p1,2)) ...
            ',' num2str(round(p2,2)) ')']);

        t=toc;
        fprintf([' (' num2str(t,2) ' sec.)']);       

        qgmdata(n).ProcessOptions = ixondata(n).ProcessOptions;
        qgmdata(n).Z = ixondata(n).Z;
        qgmdata(n).LatticeK(kk) = ixondata(n).LatticeK(kk);
        qgmdata(n).LatticePhase(kk) = ixondata(n).LatticePhase(kk);
        qgmdata(n).Params = ixondata(n).Params;
        qgmdata(n).Units = ixondata(n).Units;
        qgmdata(n).Flags = ixondata(n).Flags;
    end
    disp(' ');
end

% Caculate the average k-vector
k1Avg = k1Avg/sum(goodLatticeInds);
k2Avg = k2Avg/sum(goodLatticeInds);

% Assign the average k-vector to bad data points
if sum(goodLatticeInds)~=length(ixondata)
    warning('Bad fits for lattices found, attempting to fix using the average good ones');    
    for n=1:length(goodLatticeInds)
       if  goodLatticeInds(n)==0
           fprintf([num2str(n,'%02.f') '/' num2str(length(ixondata),'%02.f') ' ']);
             for kk=1:size(ixondata(n).Zf,3)
                % Calculate lattice phase from Fourier Transform
                ixondata(n).LatticePhase(kk) = findLatticePhase(...
                    ixondata(n).X,ixondata(n).Y,ixondata(n).Z,...
                    k1Avg,k2Avg);  
                    a1 = ixondata(n).LatticePhase(kk).a1;a2 = ixondata(n).LatticePhase(kk).a2;
                    p1 = ixondata(n).LatticePhase(kk).p1;p2 = ixondata(n).LatticePhase(kk).p2;        
                    fprintf(['a1=(' num2str(round(a1(1),2)) ...
                        ',' num2str(round(a1(2),2)) ')']);
                    fprintf(['a2=(' num2str(round(a2(1),2)) ...
                        ',' num2str(round(a1(1),2)) ')']);
                    fprintf(['phase=2*pi*(' num2str(round(p1,2)) ...
                        ',' num2str(round(p2,2)) ')']);
             end
             disp(' ');
       end
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