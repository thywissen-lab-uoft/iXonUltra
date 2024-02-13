% ixon_main_dig.m
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
%% Close all non GUI figures
% Close all figures without the GUI tag.
figs=get(groot,'Children');
disp('Closing all non GUI figures.');
for kk=1:length(figs)
   if ~isequal(figs(kk).Tag,'GUI')
        disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
        close(figs(kk)) 
   end
end
disp(' ');
%% Select image directory

    
if ~exist('dig_auto_file')
   dig_auto_file = 1; 
end

if dig_auto_file
    dialog_title='Select GUI data';       
    [filename,dig_imgdir,b]=uigetfile(fullfile(ixon_getDayDir,'*.mat'),dialog_title);
    filename = fullfile(dig_imgdir,filename);    
    if  b == 0
        disp('Canceling.');    
        return; 
    end
end
%% Load the Data
clear digdata
tic
fprintf('Loading digdata ...');
digdata = load(filename);
disp([' done (' num2str(toc,'%.2f') 's)']);

%% Initialize Options

dig_opts = struct;
dig_opts.Quality = 'auto';
dig_opts.saveDir=dig_imgdir;
strs=strsplit(imgdir,filesep);
dig_opts.FigLabel=[strs{end-1} filesep strs{end}];

%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% Choose what kind of variable to plot against (sequencer/camera)
dig_opts.varType        = 'param';          % always select 'param' for now 
dig_opts.autoXVar       = 0;                % Auto detect changing variable?
dig_opts.autoUnit       = 1;                % Auto detect unit for variable?
dig_opts.xVar           = 'ExecutionDate';  % Variable Name
dig_opts.overrideUnit   = 'V';              % If ixon_autoUnit=0, use this
dig_opts.doSave         = 1;                % Save Analysis?

%% Flags


%%
