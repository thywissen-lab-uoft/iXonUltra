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

dig_opts=bin_opts;
dig_opts.NumSigmaThresh=2.5;
digdata = bin_makeDigData2(bindata,dig_opts);

% CJF new stuff is happening.
% if exist('bin_Digitize_Source','var') && isequal(bin_Digitize_Source,'compensated')
%     digdata = bin_makeDigDataScaled(bindata,bin_opts);
% 
% else
%     bindata = ixon_digitize(bindata,dig_DigitizationThreshold);    
%     digdata = bin_makeDigData(bindata,bin_opts);
% end
    
if bin_opts.doSave 
    filename = fullfile(bin_opts.saveDir,'digdata.mat');
    disp(['Saving ' filename ' ...']);
    save(filename,'-struct','digdata');
end