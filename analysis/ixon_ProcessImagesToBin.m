function [data] = ixon_ProcessImagesToBin(data,opts)
%% Options
if nargin == 1
    opts = struct;

    % Use this if you want to use FFT to specify basis
    opts.BasisSource                = 'fft';'manual';
    opts.BasisManual                =   [0.1923 0.3244 .3;
                                        .3208 -0.1862 0.3];
    % Bin PreProcess Options
    opts.ResizeFactor               = 8;
    opts.PixelThreshold             = 150;

    % Post Binning Processing
    opts.CompensateMethod           = 'gauss';'none';'custom';

    % If radially compensating
    opts.CompensateGaussRadius      = 50;
    opts.CompensateMax              = 1.5;

    % If using a custom map
    opts.CompensateCustomMap        = 'asdfassdf.mat';

    opts.ClusterNumber              =2;
end

%% Error Checking

ProcessOptions = [data.ProcessOptions];
if (sum([ProcessOptions.doPSF])~=length(ProcessOptions));error('Cannot bin without PSF');end
if ~isfield(data,'LatticeK') && isequal(opts.BasisSource,'fft')
    error('lattice basis vectors not calculated');
end
if ~isfield(data,'LatticePhase') && isequal(opts.BasisSource,'fft')
    error('lattice basis phase not calculated');
end
if ~isfield(data,'ROI');error('No image ROI is given');end

%% Preparing
if isfield(data,'LatticeBin')
    data=rmfield(data,'LatticeBin');
end

%% Precalculate Things That Can Be Precalculated
% important improves
% speed up through ROI accumulation
% Get all ROIS
% ROI=[data(1).ROI];
% for mm = 2:length(data)
%     ROI(mm,:)=[data(mm).ROI];
% end
% 
% % If ROI is constant pre calculate things
% isOneROI= [size(unique(ROI,'rows'),1)==1];
% if isOneROI
%     X = data(1).X;
%     Y = data(1).Y;
%     ix_1 = find(X>=ROI(1),1);
%     ix_2 = find(X>=ROI(2),1);
%     iy_1 = find(Y>=ROI(3),1);
%     iy_2 = find(Y>=ROI(4),1);
%     X = X(ix_1:ix_2);
%     Y = Y(iy_1:iy_2);   
% end

%% Iterate and Bin
for kk = 1 :length(data)
    % This Stuff can be calculated ahead of time for a unique ROI
    ROI=data(kk).ROI;
    X = data(kk).X;
    Y = data(kk).Y;
    ix_1 = find(X>=ROI(1),1);
    ix_2 = find(X>=ROI(2),1);
    iy_1 = find(Y>=ROI(3),1);
    iy_2 = find(Y>=ROI(4),1);
    X = X(ix_1:ix_2);
    Y = Y(iy_1:iy_2);   
    X = linspace(X(1),X(end),length(X)*opts.ResizeFactor);
    Y = linspace(Y(1),Y(end),length(Y)*opts.ResizeFactor);    
    % Create 2xN vectors of all pixels positions
    [XX,YY]=meshgrid(X,Y);                  % matrix of X and Y
    R = [XX(:) YY(:)]'; 

    for rr = 1:size(data(kk).Z,3)        
        fprintf(['(' num2str(kk) '/' num2str(length(data)) ')' num2str(rr) ' '])
        Z = data(kk).Z(:,:,rr);

        %% Get Data
        Z = Z(iy_1:iy_2,ix_1:ix_2);   
        %% Pixel Thresholding 
        fprintf('thresh...');
        Z(Z<opts.PixelThreshold)=0;
        %% Image Resize
        % Resizing the image reduces binning issues caused from finite
        % pixel size to bin size (lattice spacing)
        fprintf('resize...');
        Z = imresize(Z,opts.ResizeFactor,'method','bilinear')/(opts.ResizeFactor^2); 
        Z = Z(:);                               %  All pixels to a list
        %% Assign Pixels To Lattice Site
        % Solve the bravais problem for every site
        fprintf('solve px2lattice ...');

        % Get the basis
        switch opts.BasisSource    
            case 'manual'
                a1 =opts.BasisManual(:,1);
                a2 =opts.BasisManual(:,2);
                p1 =opts.BasisManual(1,3);
                p2 =opts.BasisManual(2,3);
            case 'fft'
                a1 = data(kk).LatticePhase(rr).a1;
                a2 = data(kk).LatticePhase(rr).a2;                        
                p1 = data(kk).LatticePhase(rr).p1;
                p2 = data(kk).LatticePhase(rr).p2;   
        end
        A = [a1 a2];                % Basis Matrix  
        N0 = inv(A)*R;              % Pixel to Lattice Site
        p = [p1; p2];               % Phase
        P = repmat(p,[1 numel(Z)]); % Phase for all pixels
        M = round(N0-P);            % Remove phase and round to nearest site
        %% Calculate lattice vectors       
        % Solve for position of each corner in terms of the lattice basis
        na = inv(A)*[X(1); Y(1)];
        nb = inv(A)*[X(1); Y(end)];
        nc = inv(A)*[X(end); Y(1)];
        nd = inv(A)*[X(end); Y(end)];
        
        % Collect the four corners
        Nc = [na nb nc nd];
        
        % Bounds and vector for lattice basis vector 1
        n1i = floor(min(Nc(1,:)));
        n1f = ceil(max(Nc(1,:)));
        n1 = n1i:n1f;
        
        % Bounds and vector for lattice basis vector 2
        n2i = floor(min(Nc(2,:)));
        n2f = ceil(max(Nc(2,:)));
        n2 = n2i:n2f;

        %% Perform the Binning
        fprintf('binning ...');
        % Remove points outside of the desired ROI.
        ibad = logical(...
            (M(1,:)<n1i) + ...
            (M(1,:)>n1f) + ...
            (M(2,:)>n2f) + ...
            (M(2,:)<n2i));

        % Remove bad points from the binning
        M(:,ibad) = [];
        Z(ibad) = [];
        
        Zbin = zeros(n2f-n2i+1,n1f-n1i+1);
        Znum = zeros(n2f-n2i+1,n1f-n1i+1);
        for ii=1:size(M,2)    
            m1 = M(1,ii)-n1i+1;
            m2 = M(2,ii)-n2i+1;    
            Znum(m2,m1) = Znum(m2,m1)+1;
            Zbin(m2,m1) = Zbin(m2,m1) + Z(ii);  
        end
        %% Initialize Raw Binning
        ZbinRaw = Zbin;        
        %% Create Outputs
        LatticeBin                      = struct;
        LatticeBin.BinOptions           = opts;
        LatticeBin.n1                   = n1;
        LatticeBin.n2                   = n2;
        LatticeBin.ZbinRaw              = ZbinRaw;
        LatticeBin.Zbin                 = ZbinRaw;
        LatticeBin.Site2Pixel           = @(n1,n2) A*([n1;n2]+[p1;p2]);
        LatticeBin.a1                   = a1;
        LatticeBin.a2                   = a2;
        LatticeBin.p                    = [p1 p2];
        LatticeBin.LatticeSpacingPx     = mean([norm(a1) norm(a2)]);
        data(kk).LatticeBin(rr) = LatticeBin;
        disp('done');
    end
end


end

