function output = ixon_gui_bin_options

output = struct;

%% Bin Basis
% Set LatticeSource to FFT to use if you want to fit using the results from
% the FFT (recommended).
% In certain cases it is useful to use a manual specification
output.BinBasis.LatticeSource = 'manual';
output.BinBasis.BasisManual = [0.1923 0.3244 .3;
                        .3208 -0.1862 0.3];

% Use the FFT to find the basis. (Recomended)
output.BinBasis.LatticeSource = 'fft'; 

%% Binning Pre-Processing Options
% Before binning it is useful to manipulate the raw image to get better
% results
%
% Image Rescaling : Smooth issues from pixelsize~lattice spacing
% PixelThreshold  : Ignores pixels that are clearly in the noise.

output.PreBinOptions.RescaleFactor          = 8;   
output.PreBinOptions.PixelThreshold         = 150; 

%% Binning Post-Processing Options

output.PostBinOptions.RescaleRadialSigma    = 22;
output.PostBinOptions.RescaleRadialMax      = 1.2;
output.PostBinOptions.RescaleCustomFile     = 'rescale_map.mat';

%% Bin Special Function
output.StripeThreshold = [3000 7000];


end

