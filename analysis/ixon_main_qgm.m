% ixon_main_qgm.m
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
qgm_saveOpts = struct;
qgm_saveOpts.Quality = 'auto';
    
if ~exist('qgm_auto_file')
   qgm_auto_file = 1; 
end

if qgm_auto_file
    dialog_title='Select GUI data';       
    [filename,ixon_imgdir,b]=uigetfile(fullfile(ixon_getDayDir,'*.mat'),dialog_title);
    filename = fullfile(ixon_imgdir,filename);    
    if  b == 0
        disp('Canceling.');    
        return; 
    end
end

qgm_saveOpts.saveDir=ixon_imgdir;
strs=strsplit(imgdir,filesep);
qgm_saveOpts.FigLabel=[strs{end-1} filesep strs{end}];

clear qgmdata
tic
fprintf('Loading qgmdata ...');
load(filename);
disp([' done (' num2str(toc,'%.2f') 's']);

%% Match Parameters and Flags
% This makes sure that all data as has all the flags and params in each. 
% This is usually not necessary.
qgmdata = ixon_matchParamsFlags(qgmdata);

%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% Choose what kind of variable to plot against (sequencer/camera)
varType            = 'param';          % always select 'param' for now 
qgm_autoXVar       = 0;                % Auto detect changing variable?
qgm_autoUnit       = 1;                % Auto detect unit for variable?
qgm_xVar           = 'ExecutionDate';  % Variable Name
qgm_overrideUnit   = 'V';              % If ixon_autoUnit=0, use this
qgm_doSave         = 1;                % Save Analysis?

%% Flags

% Histogram for accumalted data
qgm_BinAcummulateHist                   = 1;
qgm_BinAcummulateHist_Zmax              = 6000;
qgm_BinAcummulateHist_Nbins             = 100;

% Stripe fit Data
qgm_BinStripe                           = 0;
qgm_BinStripe_LGuess                    = 25;
qgm_BinStripe_ColorThreshold            = [1000 3000];

% Digitzation
qgm_Digitize                            = 0; 
qgm_DigitizationThreshold               = 3000;

% Digization 
qgm_DigAcummulateMoments                = 0;

qgm_DigMoments                          = 0;
qgm_DigMomentsSineFit                   = 0;


%% X Variable and Units
% If auto unit and variable are chosen, search through the parameters and
% data to find which variable(s) are being changed.

% Also, sort the data by that chosen variable.

if qgm_autoXVar
    xVars = ixon_findXVars(qgmdata);
    disp([' Found ' num2str(length(xVars)) ...
        ' valid variables that are changing to plot against.']);
    disp(xVars);    
    qgm_xVar = xVars{1};    
    disp([' Setting ' qgm_xVar ' to be the x-variable']);    
    for kk=1:length(ixondata)
        disp([' (' num2str(kk) ') (' num2str(ixondata(kk).Params.(qgm_xVar)) ') ' ...
            ixondata(kk).Name]); 
    end
end

if ixon_autoUnit && isfield(ixondata(1),'Units')
    ixon_unit=ixondata(1).Units.(qgm_xVar);
else
    ixon_unit=ixon_overrideUnit;
end

if isequal(qgm_xVar,'ExecutionDate');ixon_unit='';end

% Sort the data by your given parameter
[~,inds] =  sort([P.(qgm_xVar)]);
P = [qgmdata.Params];qgmdata = qgmdata(inds);

%% Histogram

if qgm_BinAcummulateHist
    Bins = linspace(0,qgm_BinAcummulateHist_Zmax,qgm_BinAcummulateHist_Nbins);
    qgm_binnedTotalHistogram(qgmdata,Bins);
end

%% Bin Stripe
if qgm_BinStripe    
    if ~isfield(qgmdata,'LatticeBin')
        return;
    end
    clear out
    for n = 1:length(qgmdata)                
        for kk = 1:length(qgmdata(n).LatticeBin)
            n1 = qgmdata(n).LatticeBin(kk).n1;
            n2 = qgmdata(n).LatticeBin(kk).n2;
            Zb = qgmdata(n).LatticeBin(kk).Zbin;    

            opts_stripe.LGuess = qgm_BinStripe_LGuess;
            opts_stripe.FigNum=3000+10*(n-1)+kk-1;
            opts_stripe.FigNum=3000;
            opts_stripe.ColorThreshold = qgm_BinStripe_ColorThreshold;
            
            [qgmdata(n).BinStripe(kk),hF_bin_stripe] = ...
                ixon_BinStripeFit(n1,n2,Zb,opts_stripe);

            frame=getframe(hF_bin_stripe);
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);  
            
            filename = fullfile(qgm_saveOpts.saveDir,'binstripe_75.gif');
            
            if out(kk).ModDepth>=.75 && out(kk).Counts>0.5e6
                switch n
                    case 1
                        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
                    case length(qgmdata)
                        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
                    otherwise
                        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.1);
                end
            end

        end
        qgmdata(n).BinStripe = out;    
    end
end  

%% Digitization Stuff
% if ixon_doQGM_Digitize
%     ixon_main_digital;
% end