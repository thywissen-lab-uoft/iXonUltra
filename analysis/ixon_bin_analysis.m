% ixon_bin_analysis.m
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
    
if ~exist('bin_auto_file')
   bin_auto_file = 1; 
end

if bin_auto_file
    dialog_title='Select GUI data';       
    [filename,bin_imgdir,b]=uigetfile(fullfile(ixon_getDayDir,'*.mat'),dialog_title);
    filename = fullfile(bin_imgdir,filename);    
    if  b == 0
        disp('Canceling.');    
        return; 
    end
end
%% Load the Data
clear bindata
tic
fprintf('Loading bindata ...');
load(filename);
disp([' done (' num2str(toc,'%.2f') 's)']);
%% Match Parameters and Flags
% This makes sure that all data as has all the flags and params in each. 
% This is usually not necessary.
bindata = ixon_matchParamsFlags(bindata);
%% Initialize Options

bin_opts = struct;
bin_opts.Quality = 'auto';
bin_opts.saveDir=bin_imgdir;
% strs=strsplit(bin_imgdir,filesep);
% bin_opts.FigLabel=[strs{end-1} filesep strs{end}];
bin_opts.FigLabel = bindata(1).SourceDirectory;
%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% Choose what kind of variable to plot against (sequencer/camera)
bin_opts.varType        = 'param';          % always select 'param' for now 
bin_opts.autoXVar       = 0;                % Auto detect changing variable?
bin_opts.autoUnit       = 1;                % Auto detect unit for variable?
bin_opts.xVar           = 'ExecutionDate';  % Variable Name
bin_opts.overrideUnit   = 'V';              % If ixon_autoUnit=0, use this
bin_opts.doSave         = 1;                % Save Analysis?
bin_opts.DigAve          = 0;                % Digitize with average threshold? Only used with compensate thresholding

bin_opts.ControlVariable='f_offset';
% Ignore these variables when choosing auto var
autoVar_Ignore = {'f_offset','piezo_offset'};
% autoVar_Ignore = {};

%% Flags

% Histogram for accumalted data
bin_BinAcummulateHist                   = 1;
bin_FluorPerAtom                        = 1;
bin_RawVersusProcessed                  = 1;
bin_Digitize                            = 1; 

% CJF 2024/11/12 : These are now obsolete, use them at your risk
% Rescale
% bin_BinReScale                         = 1;
% Stripe fit Data
% bin_BinStripe                           = 0;
% bin_BinStripeIndex                      = 1;% 1 for trees, 2 for fallen tree
% bin_BinStripeAnimate                    = 1;
% bin_BinStripe_LGuess                    = 26.5;
% bin_BinStripe_ColorThreshold          = [3000 5000];
% bin_BinStripe_ColorThreshold            = [2500 6000];
% Digitzation
% bin_Digitize_Source                     = '
% bin_Digitize_Source                     = 'compensated';
% bin_Digitize_Source                     = 'uncompensated';



%% X Variable and Units
% If auto unit and variable are chosen, search through the parameters and
% data to find which variable(s) are being changed.

% Also, sort the data by that chosen variable.

if bin_opts.autoXVar
    xVars = ixon_findXVars(bindata);
    disp([' Found ' num2str(length(xVars)) ...
        ' valid variables that are changing to plot against.']);
    disp(xVars);    
        for kk=1:length(autoVar_Ignore)
        thisVar = autoVar_Ignore{kk};
        index = find(ismember(xVars, thisVar));
        xVars(index)=[];
    end
    
    xVar = xVars{1};    
    disp([' Setting ' xVar ' to be the x-variable']);    
    for kk=1:length(bindata)
        disp([' (' num2str(kk) ') (' num2str(bindata(kk).Params.(xVar)) ') ' ...
            bindata(kk).Name]); 
    end
    bin_opts.xVar = xVar;
end

% Sort the data by your given parameter
P = [bindata.Params];
[~,inds] =  sort([P.(bin_opts.xVar)]);
bindata = bindata(inds);

%% Standard Bin Post Processing
% After the raw bin analysis is completed from bin_initialize, there are 
% some standard things that need to be done to the binned images.
% 
% (1) Recentering : match images with same nlims (only affects displayed images)
% (2) Spatial compensation : account for fluorescence spatial inhomegeneities
% (3) intensity fluctuations : account for fluoresnece drifts

bindata = ixon_ProcessPostBin(bindata,ixon_gui_bin_options());
 
%% Histogram

if bin_BinAcummulateHist
    opts = bin_opts;
    opts.Bins       =  50;  
    opts.Nthresh    = 'auto';
    opts.saveDir    = bin_opts.saveDir;    
    opts.doAnimate  = 1;    
    opts.Source = 'Zbin';'ZbinRaw';

    % Full Cloud
    opts.ROI = 'max';
    opts.filename = 'bin_BinAnimateFull.gif';    
    hF_BinHistogramFull = bin_binnedTotalHistogram(bindata,opts);    
    if bin_opts.doSave
        ixon_saveFigure2(hF_BinHistogramFull,...
         'bin_BinHistogramFull',bin_opts);  
    end       

    % Center Cloud
    opts.ROI = [85 105 80 100];
    opts.filename = 'bin_BinAnimateCenter.gif';    
    hF_BinHistogramFull = bin_binnedTotalHistogram(bindata,opts);    
    if bin_opts.doSave
        ixon_saveFigure2(hF_BinHistogramFull,...
         'bin_BinHistogramCenter',bin_opts);  
    end    

    % Partition
    opts.NumGrid = [3 3];
    opts.Source ='ZbinRaw';
    hF_BinGridHistgoram = bin_gridHistogram(bindata,opts);    
    if bin_opts.doSave
        ixon_saveFigure2(hF_BinGridHistgoram,...
         'bin_GridHistgoram',bin_opts);  
    end    
    opts.Source ='Zbin';
end

%% Show FluorPerAtom

if bin_FluorPerAtom
    hF_fluorPeratom  = bin_showFluorPerAtom(bindata,opts);
    if bin_opts.doSave       
        ixon_saveFigure2(hF_fluorPeratom,...
            'bin_FluorPerAtom',bin_opts);  
    end 
end

%% Show Efficacy of Raw Versus Processed

if bin_RawVersusProcessed
    [hF_binCompares,output_thresh]= bin_compareHistograms(bindata,opts); 
    if bin_opts.doSave
        for kk=1:length(hF_binCompares)
            ixon_saveFigure2(hF_binCompares(kk),...
                ['bin_Compensate_' num2str(kk)],bin_opts);  
        end       
    end 
end

%% Digitization Stuff
if bin_Digitize
    ixon_dig_initialize;
end

%% Save QGM Data
if bin_opts.doSave            
    filename = fullfile(bin_opts.saveDir,'bindata.mat');
    disp(['Saving ' filename ' ...']);
    save(filename,'bindata');
end