function focus = ixon_multi_shot_focusing2(data,opts)
% Author : CF Fujiwara
%
% This code compares the focusing properties of multi shot single plane images.

disp('Two-Shot Focusing : ixon_multi_shot_focusing2');
%% Principle of Technique

% Variable which stores the piezo values, it is an array
PiezoVariableName = 'qgm_MultiPiezos';

% Indeces in the piezo array which correspond to the objective piezo value
PiezoIndeces = [3 4 5];
% PiezoIndeces = [3 4];

% Variable which stores the images (want to use unsharpenend images)
ImageSource = 'ZNoFilter';

% Indeces in the image stack to analyze
ImageIndeces = [2 3 4];
% ImageIndeces = [2 3];

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

% Z_stack_1 = zeros(size(W,1),size(W,2),length(data));
% Z_stack_2 = zeros(size(W,1),size(W,2),length(data));
% Z_stack_g = zeros(size(W,1),size(W,2),length(data));
% Zf_stack_1 = zeros(N,N,length(data));
% Zf_stack_2 = zeros(N,N,length(data));
% Zf_stack_g =  zeros(N,N,length(data));

% 
% score_1   = zeros(length(data),1);
% score_2   = zeros(length(data),1);
% score_g   = zeros(length(data),1);

% score_corr = zeros(length(data),1);

% Counts1 = zeros(length(data),1);
% Counts2 = zeros(length(data),1);

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


    % % Get the Images
    % i1      = data(kk).(ImageSource)(:,:,ImageIndex1);
    % i2      = data(kk).(ImageSource)(:,:,ImageIndex2);
    % 
    % 
    % % Blur Data
    % i1      = imgaussfilt(i1,s);
    % i2      = imgaussfilt(i2,s);
    % 
    % % Find Center of Mass
    % xc = round(sum(xx.*i1,'all')/sum(i1,'all'));
    % yc = round(sum(yy.*i2,'all')/sum(i2,'all'));
    % 
    % 
    % 
    % xR = xc + [-l:l];
    % yR = yc + [-l:l]; 
    % 
    % % Crop the Images
    % i1_crop = i1(yR,xR);
    % i2_crop = i2(yR,xR);
    % 
    % % Remove negative values
    % i1_crop(i1_crop<Nt)=0;
    % i1_crop(i1_crop<Nt)=0;
    % 
    % % Apply Window
    % i1_crop = i1_crop.*W;
    % i2_crop = i2_crop.*W;
    % 
    % N1 = sum(i1_crop,'all');
    % N2 = sum(i2_crop,'all');
    % 
    % % Normalize
    % i1_crop = i1_crop/N1;
    % i2_crop = i2_crop/N2;
    % 
    % % Find second moment on the windowed data
    % [xx2,yy2] = meshgrid([-l:l],[-l:l]);
    % xc = sum(xx2.*(i1_crop+i2_crop),'all')/sum(i1_crop+i2_crop,'all');
    % yc = sum(yy2.*(i1_crop+i2_crop),'all')/sum(i1_crop+i2_crop,'all');
    % x2 = sum(xx2.^2.*(i1_crop+i2_crop),'all')/sum(i1_crop+i2_crop,'all');
    % y2 = sum(yy2.^2.*(i1_crop+i2_crop),'all')/sum(i1_crop+i2_crop,'all');
    % 
    % % Gaussian Radius
    % xs = sqrt(x2-xc^2);
    % ys = sqrt(y2-yc^2);
    % rs = sqrt(xs.*ys);
    % 
    % % Construct fake image
    % i_gauss = exp(-(xx2-xc).^2/(2*rs^2)).*exp(-(yy2-yc).^2/(2*rs^2)).*W;
    % i_gauss = i_gauss/sum(i_gauss,'all');
    % 
    % % Calculate Peak Correlator
    % [tform,peakcorr]=imregcorr(i2_crop,i1_crop,"translation");    
    % % Rfixed  = imref2d(size(i1_crop));
    % % i2_warp = imwarp(i2_crop,tform,"OutputView",Rfixed);
    % score_corr(kk) = peakcorr;
    % 
    % % Take FFT
    % Zf_i1   = abs(fftshift(fft2(i1_crop)));
    % Zf_i2   = abs(fftshift(fft2(i2_crop)));
    % Zf_ig   = abs(fftshift(fft2(i_gauss)));
    % 
    % % Radial Sum of FFT (for plotting purposes)
    % [r,Zf_i1_r,~,~] = radial_profile(Zf_i1,nBin);
    % [r,Zf_i2_r,~,~] = radial_profile(Zf_i2,nBin);
    % [r,Zf_ig_r,~,~] = radial_profile(Zf_ig,nBin);
    % fr      = r/N;    % Radial Frequency Vector 1/px
    % 
    % if kk==1
    %     Zfr_stack_1 = zeros(length(r),length(data));
    %     Zfr_stack_2 = zeros(length(r),length(data));
    %     Zfr_stack_g = zeros(length(r),length(data));
    % end
    % 
    % % Cartesian Frequency Vector 1/px
    % f1      = linspace(-0.5,0.5,size(Zf_i1,1));
    % f2      = linspace(-0.5,0.5,size(Zf_i1,2));
    % 
    % % Compute Momentum Score
    % [ffx,ffy]=meshgrid(f2,f1);
    % ffr     = sqrt(ffx.^2+ffy.^2);
    % s1      = sum(Zf_i1.*ffr,'all');
    % s2      = sum(Zf_i2.*ffr,'all');
    % sg      = sum(Zf_ig.*ffr,'all');
    % 
    % % Make Outputs
    % Z_stack_1(:,:,kk) = i1_crop*N1;
    % Z_stack_2(:,:,kk) = i2_crop*N2;
    % Z_stack_g(:,:,kk) = i_gauss*(N1+N2)*0.5;
    % 
    % Zfr_stack_1(:,kk) = Zf_i1_r/nBin;
    % Zfr_stack_2(:,kk) = Zf_i2_r/nBin;
    % Zfr_stack_g(:,kk) = Zf_ig_r/nBin;
    % 
    % Zf_stack_1(:,:,kk) = Zf_i1;
    % Zf_stack_2(:,:,kk) = Zf_i2;
    % Zf_stack_g(:,:,kk) = Zf_ig;
    % 
    % Counts1(kk) = N1;
    % Counts2(kk) = N2;
    % 
    % score_1(kk) = s1;
    % score_2(kk) = s2;
    % score_g(kk) = sg;

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


% 
% focus.dSdV          = ((score_1-score_2)./score_g)./(piezos(:,1)-piezos(:,2));
% 
% focus.Image1            = Z_stack_1;
% focus.Image2            = Z_stack_2;
% % focus.ImageG            = Z_stack_g;
% % 
% focus.Image1_FFT        = Zf_stack_1;
% focus.Image2_FFT        = Zf_stack_2;
% % focus.ImageG_FFT        = Zf_stack_g;
% 
% focus.Image1_FFT_radial = Zfr_stack_1;
% focus.Image2_FFT_radial = Zfr_stack_2;
% focus.ImageG_FFT_radial = Zfr_stack_g;
% 
% focus.BoxCount1         = Counts1;
% focus.BoxCount2         = Counts2;
% 
% focus.Piezo1            = piezos(:,1);
% focus.Piezo2            = piezos(:,2);
% focus.Score1            = score_1;
% focus.Score2            = score_2;
% focus.ScoreGauss        = score_g;
% focus.Correlator        = score_corr;


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