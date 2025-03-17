function [data,focus_data,hF,hF_Summary] = ixon_multi_shot_focusing2(data,opts)
% Author : CF Fujiwara
%
% This code compares the focusing properties of subsequent images
% 
% It assume that the image stack comprises of three images where image 2
% and image 3 are taken at different focus positions of the microscope.  We
% do not use the 1st image as we allow for atoms leave to exit the trap
% just in case they are not properly cooled. While the analysis compares
% the normalized count profile between the images, it works best if the
% light source (ie. atomic distribution) is identical from shot to shot.

focus_data=struct;
hF=[];
hF_Summary=[];


% This variable stores the piezo values, it is an array
PiezoVariableName = 'qgm_MultiPiezos';
PiezoIndex1 = 3;
PiezoIndex2 = 4;


ImageSource = 'ZNoFilter';
ImageIndex1 = 2;
ImageIndex2 = 3;

%% Get Piezo Values
% In variable exposure mode, there are a total of eight exposures.
% [Wipe Image1 Image2 Image3 Wipe Image4 Image5 Image6]
%
% As the code is updated, this function may need to change.

piezos = zeros(length(data),2);
P=[data.Params];
for kk=1:length(data)
    vals=[P(kk).(PiezoVariableName)]; 
    piezos(kk,1)=vals(PiezoIndex1);
    piezos(kk,2)=vals(PiezoIndex2);
end

%% Default Settings

if nargin==1; opts=struct;end

% Default ROI
% When choosing an ROI, it is important to NOT allow for clipping induced
% by the image plane iris or any other "hard" features.
if ~isfield(opts,'ROI')
    opts.ROI = [150 250 150 250];
end

% Images to analyze
if ~isfield(opts,'ImageIndeces')
    opts.ImageIndeces = [2 3];
end

%%
ROI = opts.ROI;

for kk=1:length(data)
    % Get the Images
    i1      = data(kk).(ImageSource)(:,:,ImageIndex1);
    i2      = data(kk).(ImageSource)(:,:,ImageIndex2);

    % Crop the Images
    i1_crop = i1(ROI(3):ROI(4),ROI(1):ROI(2));
    i2_crop = i2(ROI(3):ROI(4),ROI(1):ROI(2));

    % Calculate Peak Correlator
    [tform,peakcorr]=imregcorr(i2_crop,i1_crop,"translation");    
    Rfixed  = imref2d(size(i1_crop));
    i2_warp = imwarp(i2_crop,tform,"OutputView",Rfixed);

    % Take fft
    Zf_i1   = abs(fftshift(fft2(i1_crop)));
    Zf_i2   = abs(fftshift(fft2(i2_crop)));
    f1      = linspace(-1,1,size(Zf_i1,1));
    f2      = linspace(-1,1,size(Zf_i1,2));

    [ffx,ffy]=meshgrid(f2,f1);
    ffr=sqrt(ffx.^2+ffy.^2);
    s1(kk)=sum(Zf_i1.*ffr,'all');
    s2(kk)=sum(Zf_i2.*ffr,'all');


    x = ixondata(n).X;
    y = ixondata(n).Y;
end

end



