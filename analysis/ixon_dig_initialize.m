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

%%

qgm_DigitizationThreshold               = 3500;
qgmdata = ixon_digitize(qgmdata,qgm_DigitizationThreshold);    
digdata = qgm_makeDigData(qgmdata,qgm_opts);
    
if qgm_opts.doSave 
    filename = [qgm_opts.saveDir 'digdata.mat'];
    save(filename,'-struct','digdata');
end