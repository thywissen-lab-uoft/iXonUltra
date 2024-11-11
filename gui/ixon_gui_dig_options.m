function opts = ixon_gui_dig_options

opts=struct;

%% Digitization Method Options

% Use fitted probability density function [pdf] to find the threshold
opts.DigitizeMethod     = 'pdf_kmeans_threshold'; % recommended (uses kmeans+pdf)
opts.pdf_Sigma          = 2.5;

% % Use kmeans clustering to find threshold
% opts.DigitzeMethod      = 'kmeans';
% 
% % Use a manual threshold
% opts.DigitizeMethod     = 'manual';
% opts.ManualThreshold    = 2000;

% For future useful to not just do simple thresholding by somekind of
% algorithm to calculate
% opts.MaximumLikeliHoodAlgorithm = false;

%% Digital Trap Calibrations
opts.doHubbardAnalysis = true;


% Trap Frequency
% The overall harmonic trap frequency
opts.TrapOmega      = 2*pi*65;

% Lattice Depth V0 [Er]
% In analyzing the images, the lattice depth tells you the 
opts.LatticeDepth   = 2.5; % For Manual specification of depth
opts.LatticeDepth   = 'lattice_depth_var_name';    % Pull depth from these params

% Tunneling t [Hz]
% The single particle tunneling element.
opts.Tunneling = [];    % keep empty if auto-calculate from V0 
opts.Tunneling = 563;

% Magnetic Field
% The magnetic field is used to calculate the s-wave scattering length.
opts.MagneticField  = 'varname';    % 
opts.MagneticField  = 201.5;        % Specify magnetic field manually

% Wannier Source
opts.WannierOverlap = [];           % keep empty if auto-calculate from V0

% Interaction Strength U [Hz]
% The hubbard can calculated from V0 and a_s.  In the perturbative limit
% the hubbard U is given by the WannierOverlap function.
opts.U = [];    % Keep this empty for automatic calculation



%% Digital Analysis Options

% Gauss Fitting Options
% Fit the radial distribution to a radial gaussian.  To account for pauli
% and interactino effects, only fit the wings where the density lower.
opts.GaussFitDensityMax = [];       % empty if you don't want gauss fit
opts.GaussFitDensityMax =[0.05 1];  % Maximum density to fit to gaussian


end

