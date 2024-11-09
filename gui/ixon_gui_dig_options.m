function opts = ixon_gui_dig_options

opts=struct;

%% Digitization Method Options

% Use fitted probability density function [pdf] to find the threshold
opts.DigitizeMethod     = 'pdf_kmeans_threshold'; % recommended (uses kmeans+pdf)
opts.pdf_Sigma          = 2.5;

% Use kmeans clustering to find threshold
opts.DigitzeMethod      = 'kmeans';

% Use a manual threshold
opts.DigitizeMethod     = 'manual';
opts.ManualThreshold    = 2000;

% For future useful to not just do simple thresholding by somekind of
% algorithm to calculate
opts.MaximumLikeliHoodAlgorithm = false;

%% Analysis Options

% Trap Frequency
opts.TrapOmega      = 2*pi*[65 65];

% Lattice Depth/Tunneling
opts.LatticeDepth   = [2.5 2.5 2.5]; % For Manual specification of depth
opts.LatticeDepth   = {'l','1','1'};    % Pull depth from these params

% Tunneling Source
% How to convert lattice depth in tunneling

% Wannier Source
% How to calculate the wannier interaction overlap

% Magnetic Field (interactions)
opts.MagneticField  = 201.5;
opts.MagneticField  = 'varname';        % Pull field from this param

end

