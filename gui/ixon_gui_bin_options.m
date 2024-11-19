function opts = ixon_gui_bin_options

%% Options for binning
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

