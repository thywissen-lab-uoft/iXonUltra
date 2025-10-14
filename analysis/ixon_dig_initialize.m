% ixon_main_dig_initialize.m
%
% Author : CF Fujiwara
%
% This script is the primary analysis code for ixon images which are single
% plane.

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));    

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath));

a = fileparts(curpath);
addpath(a);addpath(genpath(a));

%% Description
% Once all of the binning has been completed, the next stage of the
% analysis comes to the digitization. Likewise with the binning protocol,
% there is quite a large range of complicated protocols that can be
% applied.  In most cases, you can just run the code, but I (CJF) highly
% recommends that you fully understand the pros/cons of the various
% protocols.
%
% After binning, the following information should be known :
%   - All fluoresnce counts are binned into lattice sites
%   - Images with few atoms have lattices wavevectors chosen based upon
%   other images in the data set
%   - All fluorescence counts are scaled based on there distance from the
%   center of the cloud. This is mean to account for inhomogeneities in the
%   fluoresnce pattern.
%   - a statistical analysis of the binned counts is done.

%% Initialize Digdata

dig_opts=bin_opts;
dig_opts.NormalizedCenter = 1.07;       % Center point of normalized distribution (close to unity)
dig_opts.NormalizedSigma  = 0.24;       % Gaussian radius of normalized distribution
dig_opts.NormalizedThreshold = ...      % Simple thresholding
    dig_opts.NormalizedCenter - 1*dig_opts.NormalizedSigma;

digdata = bin_makeDigData2(bindata,dig_opts);


if bin_opts.doSave 
    filename = fullfile(bin_opts.saveDir,'digdata.mat');
    disp(['Saving ' filename ' ...']);
    save(filename,'-struct','digdata');
end