function opts = ixon_gui_K_options

opts=struct;

%% Lattice Basis and Phase Options

% CJF: Currently not used. Would need to modify code to do this, but this
% is a sensible change

%% Multi-Shot Focusing Options

opts.KFocus = struct;
opts.KFocus.kmag = 0.3727;
opts.KFocus.kdelta = 0.01;
opts.KFocus.ControlVariable = 'qgm_MultiPiezos';
opts.KFocus.doDebug = 1;
opts.KFocus.ROI = [1 512 1 512];

end

