function focus = ixon_MultiShotFocus(data,opts)
% Author : CF Fujiwara
%
% This code compares the focusing properties of multi shot single plane images.

if nargin == 1
    opts = struct;
end
disp('Two-Shot Focusing : ixon_multi_shot_focusing2');
%% Default Settings

if ~isfield(opts,'PiezoVariableName')
    opts.PiezoVariableName = 'qgm_MultiPiezos';
end

if ~isfield(opts,'PiezoIndeces')
    opts.PiezoIndeces = [3 4 5];
end

if ~isfield(opts,'ImageIndeces')
    opts.ImageIndeces = [2 3 4];
end

if ~isfield(opts,'ImageSource')
    opts.ImageSource = 'ZNoFilter';
end
%% Store Values

% Variable which stores the piezo values, it is an array
PiezoVariableName = opts.PiezoVariableName;
PiezoIndeces = opts.PiezoIndeces;
ImageSource = opts.ImageSource;
ImageIndeces = opts.ImageIndeces;

disp(['     PiezoIndeces   : ' num2str(PiezoIndeces)]);
disp(['     ImageIndeces   : ' num2str(ImageIndeces)]);

%% Get Piezo Values

piezos = zeros(length(data),length(PiezoIndeces));
P=[data.Params];
for kk=1:length(data)
    vals=[P(kk).(PiezoVariableName)]; 
    for jj=1:length(PiezoIndeces)
        piezos(kk,jj)=vals(PiezoIndeces(jj));
    end
end

D=[P.ExecutionDate]';
%% Default Settings

if nargin==1; opts=struct;end

% Images to analyze
if ~isfield(opts,'ImageIndeces')
    opts.ImageIndeces = [2 3];
end

% Blur radius
s = 1.5;

% Window Radius
l = 150;

% Threshold
Nt = 0;

% Radial FFT bin
nBin = 6;
Nfft = 2*l+1;

disp(['     Blur Radius    : ' num2str(s)]);
disp(['     Hanning Radius : ' num2str(l)]);
%% Create Window

w1 = hann(2*l+1);
w2 = hann(2*l+1);
W = w1 * w2';

%% Prepare output Data

% Position Domain Image Stack [row col piezo run#]
Z_stack  = zeros(size(W,1),size(W,2),length(data),length(PiezoIndeces));

% Frequency Domain Image Stack [row col piezo run#]
Zf_stack = zeros(size(W,1),size(W,2),length(data),length(PiezoIndeces));

% Box Counts
counts   = zeros(length(data),length(PiezoIndeces));

% Focus Scores
scores      = zeros(length(data),length(PiezoIndeces));

% Correlation Scores
score_corr  = zeros(length(data),length(PiezoIndeces)-1);


%%

for kk=1:length(data)
    fprintf([num2str(kk) ' of ' num2str(length(data)) ' ...']);
    tic;

    % Get Position Vector
    x       = data(kk).X;
    y       = data(kk).Y;
    [xx,yy]=meshgrid(x,y);

    for jj=1:length(PiezoIndeces)
        img = data(kk).(ImageSource)(:,:,ImageIndeces(jj));
        img = imgaussfilt(img,s);
    
        % x center of mass
        xc = round(sum(xx.*img,'all')/sum(img,'all'));
        xc = max([xc l+1]);
        xc = min([xc 512-l-1]);

        % Y center of mass
        yc = round(sum(yy.*img,'all')/sum(img,'all'));
        yc = max([yc l+1]);
        yc = min([yc 512-l-1]);

        % Sub-Vector
        xR = xc + [-l:l];
        yR = yc + [-l:l]; 

        img_crop    = img(yR,xR);           % Crop Image        
        img_crop(img_crop<0) = 0;           % Remove Negative Value
        img_crop    = img_crop.*W;          % Apply window
        N           = sum(img_crop,'all');  % Normalize        
        img_crop    = img_crop/N;           % Normalize
        img_fft     = ...                   % Take FFT
            abs(fftshift(fft2(img_crop))); 
        [r,img_fft_r,~,~] = ...               % Radial FFT
            radial_profile(img_fft,nBin);
        fr      = r/Nfft;    % Radial Frequency Vector 1/px

        if kk==1 && jj==1
            Zfr_stack = zeros(length(r),length(data),length(PiezoIndeces));
            % Cartesian Frequency Vector 1/px
            f1      = linspace(-0.5,0.5,size(W,1));
            f2      = linspace(-0.5,0.5,size(W,2));
            [ffx,ffy]=meshgrid(f2,f1);
            ffr     = sqrt(ffx.^2+ffy.^2);
        end

        counts(kk,jj) = N;
        scores(kk,jj) = sum(img_fft.*ffr,'all');
      
        Zfr_stack(:,kk,jj) = img_fft_r(:);
        Z_stack(:,:,kk,jj) = img_crop;
        Zf_stack(:,:,kk,jj) = img_fft;
        
    end

    % Calculate Peak Correlator
    % [tform,peakcorr]=imregcorr(i2_crop,i1_crop,"translation");    
    % % Rfixed  = imref2d(size(i1_crop));
    % % i2_warp = imwarp(i2_crop,tform,"OutputView",Rfixed);
    % score_corr(kk) = peakcorr;   

    t2=toc;
    disp([' done (' num2str(round(1e3*t2)) ' ms)'])
end

%% Construct Output

focus = struct;
focus.Scores            = scores;
focus.Counts            = counts;
focus.RadialFrequency   = fr;
focus.RadialFFT         = Zfr_stack;
focus.Images_Pos        = Z_stack;
focus.Images_Freq       = Zf_stack;
focus.Piezos            = piezos;
focus.ExecutionDate     = D;
focus.Params            = P;
focus.PiezoVariableName = 'qgm_MultiPiezos';
focus.PiezoIndeces      = PiezoIndeces;
focus.ImageSource       = ImageSource;
focus.ImageIndeces      = ImageIndeces;
focus.BlurRadius        = s;
focus.WindowRadius      = l;
focus.FreqnBin          = nBin;

end

function [Tics,Average,dev,n]=radial_profile(data,radial_step)
    %main axii cpecified:
    x=(1:size(data,2))-size(data,2)/2;
    y=(1:size(data,1))-size(data,1)/2;
    % coordinate grid:
    [X,Y]=meshgrid(x,y);
    % creating circular layers
    Z_integer=round(abs(X+1i*Y)/radial_step)+1;
    % % illustrating the principle:
    % % figure;imagesc(Z_integer.*data)
    % very fast MatLab calculations:
    Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
    Average=accumarray(Z_integer(:),data(:),[],@sum);
    
    
    % dev=accumarray(Z_integer(:),data(:),[],@std);
    dev=accumarray(Z_integer(:),data(:),[],@(x) std(x,1));
    n= accumarray(Z_integer(:),data(:),[],@(x) numel(x));
end