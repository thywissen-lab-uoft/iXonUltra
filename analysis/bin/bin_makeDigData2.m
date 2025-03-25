function digdata = bin_makeDigData2(bindata,opts)
% Author : CJ Fujiwara
% Last Edited : 2025/03/25
%% Settings
if nargin ==1
    opts=struct;
end

if ~isfield(opts,'xVar')
    opts.xVar='ExecutionDate';
end

if ~isfield(opts,'NumSigmaThresh')
    opts.NumSigmaThresh=2.5;
end

if ~isfield(opts,'NormThresh')
    opts.NormThresh=[];
end

if ~isfield(opts,'CountMax')
    opts.CountMax = 35e3;
end

%% Grab Basic Data Information

P = [bindata.Params];           % Params structure
F = [bindata.Flags];            % Flags structure
U = [bindata.Units];            % Units struture
B = [bindata.LatticeBin];       % Bin Data
PH =[bindata.LatticePhase];     % Phase Data

%% Initialize Data Structures
Zdig = zeros(...
    size(bindata(1).LatticeBin(1).Zbin,1),...
    size(bindata(1).LatticeBin(1).Zbin,2),...
    length(bindata),...
    length(bindata(1).LatticeBin));
% The digitized data is MxNxLxR. 
% MxN is size of image
% L is length of runs
% R is number of data images in each run

% Lattice vectors
n1 = bindata(1).LatticeBin.n1;
n2 = bindata(1).LatticeBin.n2;

% Convert lattice vectors in 2D array and then go back to 1D
[nn1,nn2]=meshgrid(n1,n2);
N1_all = nn1(:);
N2_all = nn2(:);
    
%% Acquire all thresholds
for kk=1:length(bindata)
    for rr=1:length(bindata(kk).LatticeBin)
        if ~isempty(bindata(kk).LatticeBin(rr).PDF1_Center)
            count_center(kk,rr) = bindata(kk).LatticeBin(rr).PDF1_Center;
            count_sigma(kk,rr) = bindata(kk).LatticeBin(rr).PDF1_Radius;
        else
            % Its empty if there are not enough atoms to measure it
            count_center(kk,rr) = NaN;
            count_sigma(kk,rr) = NaN;    
        end            
    end
end

% Average over each image
count_center=mean(count_center,2);
count_sigma=mean(count_sigma,2);

% Find bad indeces
bad_inds = isnan(count_center);

% Remove bad thresholds
count_center_good = count_center;
count_center_good(bad_inds)=[];
count_sigma_good = count_sigma;
count_sigma_good(bad_inds)=[];

% Bad Threshods becomes the average value
count_center(bad_inds) = mean(count_center_good,'all');
count_sigma(bad_inds) = mean(count_sigma_good,'all');

count_threshold = opts.NormalizedThreshold*count_center;


    %% Do the Digitization

    % Initialize Stuff, This seems excessive
    Site2PixelsFuncs={};
    Siteatom={};
    Ratom = {};
    Natoms=[];
    lattice_spacing_px=[];
    Counts={};


    for kk=1:length(bindata)
        Site2PixelsFuncs{kk}=bindata(kk).LatticeBin(1).Site2Pixel;
        a1=norm(bindata(kk).LatticeBin(rr).a1);
        a2=norm(bindata(kk).LatticeBin(rr).a1);
        lattice_spacing_px(kk) = mean([a1 a2]);
        for rr=1:length(bindata(kk).LatticeBin)
            Zbin_rr    = bindata(kk).LatticeBin(rr).Zbin;   % Binned Data            
            Zdig_rr    = [Zbin_rr>=count_threshold(kk)];    % Simple Thresholding     
            Zdig(:,:,kk,rr) = Zdig_rr;                      % Add digital matrix to data            
            
            % List of all sites that are occupied
            N1 = N1_all;                            % All sites N1
            N2 = N2_all;                            % All sites N2
            N1(~Zdig_rr(:)) = [];                   % Remove Empty Sites from N1
            N2(~Zdig_rr(:)) = [];                   % Remove Empty Sites from N2
            N1=N1';N2= N2';                         % Turn into 1xL

            % Vector with all occupied sites
            Siteatom{kk,rr} = [N1;N2];      
            % Vector with all occupied pixel positions
            Ratom{kk,rr}    = bindata(kk).LatticeBin(rr).Site2Pixel(N1,N2);
            % Number of atoms in an image
            Natoms(kk,rr)   = sum(Zdig_rr,'all');
            % Fluoresence counts attributed to atoms
            Counts{kk,rr}   = Zbin_rr(Zdig_rr);  
        end
    end        
    
    
   %% Output    
    digdata                     = struct;    
    digdata.SourceDirectory     = unique({bindata.SourceDirectory});
    digdata.FileNames           = {bindata.Name}';
    digdata.Threshold           = count_threshold;
    digdata.ThresholdType       = 'counts';
    digdata.xVar                = opts.xVar;
    digdata.X                   = [P.(opts.xVar)];
    digdata.Params              = P;
    digdata.Units               = U;
    digdata.Flags               = F;      
    
    % Actual Digital Data
    digdata.n1                  = n1;           % site vector
    digdata.n2                  = n2;           % site vector
    digdata.Zdig                = logical(Zdig);% digitized image
    digdata.Counts              = Counts;       % counts for each atom
    digdata.Ratom               = Ratom;        % pixel position of each atom
    digdata.Siteatom            = Siteatom;     % site index of each atom
    digdata.Natoms              = Natoms;       % number of atoms
    digdata.Lattice_px          = mean(lattice_spacing_px); % lattice spacing
    digdata.Lattice_um          = 1.054/2;      % lattice spacing

    % Lattice Spacing

    % Basis Information
    digdata.a1                  = [B.a1];       % Lattice Vector 1
    digdata.a2                  = [B.a2];       % Lattice Vector 2    
    digdata.p1                  = [PH.p1];      % Lattice Phase 1
    digdata.p2                  = [PH.p2];      % Lattice Phase 2
    
    
    [digdata] = dig_basic(digdata);