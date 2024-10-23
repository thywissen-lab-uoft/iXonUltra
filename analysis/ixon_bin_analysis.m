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

bin_opts.ControlVariable='f_offset';

%% Flags

% Recenter all binned data to have same limits
bin_BinRecenter                         = 1;

% Histogram for accumalted data
bin_BinAcummulateHist                   = 1;
bin_BinAcummulateHist_Zmax              = 30000;
bin_BinAcummulateHist_Nbins             = 100;

% Rescale
bin_BinReScale                          = 1;


% Stripe fit Data
bin_BinStripe                           = 1;
bin_BinStripeAnimate                    = 1;
bin_BinStripe_LGuess                    = 26.5;
bin_BinStripe_ColorThreshold            = [3000 5000];

% Digitzation
bin_Digitize                            = 0; 
dig_DigitizationThreshold               = 2000;

bin_Digitize_Source                     = 'compensated';
%bin_Digitize_Source                     = 'uncompensated';

%% X Variable and Units
% If auto unit and variable are chosen, search through the parameters and
% data to find which variable(s) are being changed.

% Also, sort the data by that chosen variable.

if bin_opts.autoXVar
    xVars = ixon_findXVars(bindata);
    disp([' Found ' num2str(length(xVars)) ...
        ' valid variables that are changing to plot against.']);
    disp(xVars);    
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

%% Recenter

if bin_BinRecenter
    bindata = bin_recenter(bindata);
end

 
%% Histogram

if bin_BinAcummulateHist
    opts = bin_opts;
    opts.Bins =  linspace(0,bin_BinAcummulateHist_Zmax,...
    bin_BinAcummulateHist_Nbins);  
    opts.Nthresh =dig_DigitizationThreshold;

    opts.saveDir = bin_opts.saveDir;
    
    opts.doAnimate = 1;
    
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
    hF_BinGridHistgoram = bin_gridHistogram(bindata,opts);    
    if bin_opts.doSave
        ixon_saveFigure2(hF_BinGridHistgoram,...
         'bin_GridHistgoram',bin_opts);  
    end
end

%%

if bin_BinReScale
    [bindata,hF_BinCompensate] = bin_rescale(bindata,opts);    
    if bin_opts.doSave
        ixon_saveFigure2(hF_BinCompensate,...
         'bin_Compensate',bin_opts);  
    end
    
end
%% Bin Stripe
if bin_BinStripe    
    if ~isfield(bindata,'LatticeBin')
        return;
    end
    clear out
    if isfield(bindata,'BinStripe')
    bindata=rmfield(bindata,'BinStripe');
    end
    for n = 1:length(bindata)    
        fprintf(['(' num2str(n) '/' num2str(numel(bindata))...
            ') lattice stripe fit']);
        tic
        for kk = 1:length(bindata(n).LatticeBin)
            fprintf(['...' num2str(kk)]);
            n1 = bindata(n).LatticeBin(kk).n1;
            n2 = bindata(n).LatticeBin(kk).n2;
            Zb = bindata(n).LatticeBin(kk).Zbin(:,:,1);    

            opts_stripe.LGuess = bin_BinStripe_LGuess;
            opts_stripe.FigNum=3000+10*(n-1)+kk-1;
            opts_stripe.FigNum=3000;
            opts_stripe.Threshold = bin_BinStripe_ColorThreshold;
            
            bindata(n).BinStripe(kk) = ...
                bin_StripeFit(n1,n2,Zb,opts_stripe);
        end
        disp([' done (' num2str(toc,'%.2f') 's)']); 
    end

end  
%% Bin Stripe Summary
if bin_BinStripe    
    bin_opts.nCenter = 110;
    hF_StripeSummary = bin_showStripeBinSummary(bindata,bin_opts.xVar,bin_opts);    
    if bin_opts.doSave
        ixon_saveFigure2(hF_StripeSummary,...
         'bin_StripeSummary',opts);     
    end
end

%% Bin Stripe Animation
if bin_BinStripe && bin_BinStripeAnimate
    opts = bin_opts;
    opts.Threshold = bin_BinStripe_ColorThreshold;
    opts.filename = 'bin_BinStripeAnimation.gif';
    opts.SumIndex = 1; % 1 for trees % 2 for fallen trees

    bin_showStripeBin(bindata,bin_opts.xVar,opts);
end
%% Digitization Stuff
if bin_Digitize
    ixon_dig_initialize;
end