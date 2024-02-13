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

%% Initialize Digdata

dig_DigitizationThreshold               = 3500;
bindata = ixon_digitize(bindata,dig_DigitizationThreshold);    
digdata = bin_makeDigData(bindata,bin_opts);
    
if bin_opts.doSave 
    filename = [bin_opts.saveDir 'digdata.mat'];
    disp(['Saving ' filename ' ...']);
    save(filename,'-struct','digdata');
end