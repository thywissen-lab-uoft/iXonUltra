function opts = ixon_gui_dig_options

opts=struct;

%% Digitization Method Options

% Use fitted probability density function [pdf] to find the threshold
opts.DigitizeMethod     = 'pdf_kmeans'; % recommended (uses kmeans+pdf)
opts.pdf_Sigma          = 2.5;

% Use kmeans clustering to find threshold
opts.DigitzeMethod      = 'kmeans';

% Use a manual threshold
opts.DigitizeMethod     = 'manual';
opts.ManualThreshold    = 2000;

%% Analysis Options

opts.TrapOmega      = 2*pi*[65 65];
opts.Tunneling      = [563 563 563];
opts.MagneticField  = 201.5;

end

