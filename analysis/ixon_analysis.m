% ixonAnalysis.m
%
% Author : CF Fujiwara
%
% This script is the primary analysis file for the iXon MATLAB code. It
% calls and plots all other analyses.

% Display this filename
disp(repmat('-',1,60));    
disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));    
disp(repmat('-',1,60));    

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath));
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

%% Select global settings
% This section of code sets some global settings for the analysis which the
% imaging GUI is unaware of.  

disp('Choosing global settings for analysis...');

global ixon_imgdir

lambda=770E-9;

% Load pertinent physical constants
amu=1.660539E-27; 
m=40*amu;

%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% Choose what kind of variable to plot against (sequencer/camera)
varType='param'; % always select 'param' for now 


ixon_xVar='ExecutionDate';
unit='s';

% Flag whether to save the output figures or not (code is faster if not
% saving)
ixon_doSave=1;

%% Select image directory
% Choose the directory where the images to analyze are stored
disp([datestr(now,13) ' Choose an image analysis folder...']);
dialog_title='Choose the root dire ctory of the images';
ixon_imgdir=uigetdir(ixon_getImageDir(datevec(now)),dialog_title);
if isequal(ixon_imgdir,0)
    disp('Canceling.');
    return 
end

%% Load the data
clear ixondata
disp(['Loading data from ' ixon_imgdir]);
files=dir([ixon_imgdir filesep '*.mat']);
files={files.name};


for kk=1:length(files)
    str=fullfile(ixon_imgdir,files{kk});
    [a,b,c]=fileparts(str);      
    disp(['     (' num2str(kk) ')' files{kk}]);    
    data=load(str);     
    data=data.data;  

    % Display image properties
    try
        disp(['     Image Name     : ' data.Name]);
        disp(['     Execution Time : ' datestr(data.Date)]);
        disp(['     ' ixon_xVar ' : ' num2str(data.Params.(ixon_xVar))]);
        disp(' ');
    end    
    
    if isequal(ixon_xVar,'ExecutionDate')
        data.Params.(ixon_xVar)=datenum(data.Params.(ixon_xVar))*24*60*60;
    end
    
    ixondata(kk)=data;    
end
disp(' ');

if isequal(ixon_xVar,'ExecutionDate')
   p=[ixondata.Params] ;
   tmin=min([p.ExecutionDate]);
   for kk=1:length(ixondata)
      ixondata(kk).Params.ExecutionDate= ...
          ixondata(kk).Params.ExecutionDate-tmin;
   end     
end

%% Sort the data
% Sort the data by your given parameter
clear x
disp(['Sorting ixondata by the given ''' ixon_xVar '''']);

% Get the object that contains the variable
switch varType
    case 'param'
        varList=[ixondata.Params];
    case 'acq'
        varList=[ixondata.AcquisitionInformation];        
    otherwise
        error('uhh you chose the wrong thing to plot');
end
       
% Make sure that all ixondata have the parameter and record it
allGood=1;
for kk=1:length(varList)
    if isfield(varList(kk),ixon_xVar)
        x(kk)=varList(kk).(ixon_xVar);
    else
        allGood=0;
        warning(['ixondata(' num2str(kk) ') has no ''' ixon_xVar '''']);
    end
end

% Sort if all the data has that parameter
if allGood
    [~, inds]=sort(x);
    ixondata=ixondata(inds);
end

% Get the object that contains the variable
switch varType
    case 'param'
        varList=[ixondata.Params];
    case 'acq'
        varList=[ixondata.AcquisitionInformation];        
    otherwise
        error('uhh you chose the wrong thing to plot');
end


%% Analysis ROI
% Analysis ROI is an Nx4 matrix of [X1 X2 Y1 Y2] which specifies a region
% to analyze. Each new row in the matrix indicates a separate ROI to
% perform analysis on.
%
% While in principle different images can have different analysis ROIs,
% this is currently disabled because it creates code issues at the moment.

% Full ROI
ixonROI = [1 512 1 512];   


[ixondata.ROI]=deal(ixonROI);

%% Image Processing

% Electronic bias offset
doSubBias=1;
offset=200; % (Doesnt affect calculations)

% Ixon mask
doApplyMask=1;
maskname=fullfile('ixon_mask.mat');
ixon_mask=load(maskname);
ixon_mask=ixon_mask.BW;

% Gauss Filter
doGaussFilter=1;
filter_radius=0.25;  % Gaussian filter radius


for kk=1:length(ixondata)
    imgs=ixondata(kk).RawImages;    
    
    for jj=1:size(ixondata(kk).RawImages,3)   
        if doSubBias
           imgs(:,:,jj)=imgs(:,:,jj)-offset; 
        end
        
        if doApplyMask
           imgs(:,:,jj)=imgs(:,:,jj).*ixon_mask;
        end
        
        if doGaussFilter
           imgs(:,:,jj)=imgaussfilt(imgs(:,:,jj),filter_radius);
        end
    end
    
    % For now we assume only two images
    Z=imgs(:,:,2)-imgs(:,:,1);
    ixondata(kk).Z=Z;
end



%% Basic Raw Image Analysis

doRawImageAnalysis=0;
if doRawImageAnalysis   

    % Do basic analysis on raw counts
    ixondata=ixon_computeRawCounts(ixondata);

    % Plot histogram of raw counts
    hist_opts=struct;
    hist_opts.Outliers=[10 50]; % Histogram wont plot outliers of this many low/high
    hist_opts.GlobalLimits=1;   % Maintain historgram x limits
    hist_opts.BinWidth=10;       % Histogram bin width
    hist_opts.ImageNumber=1;    % Which image to histogram (overwritten)
    hist_opts.YScale='Log';     % Histogram y scale
    % hist_opts.YScale='Linear';

    for kk=1:size(ixondata(1).RawImages,3)
        hist_opts.ImageNumber=kk;
        hF_ixon_rawhist=ixon_showRawCountHistogram(ixondata,ixon_xVar,hist_opts);
        if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_rawhist,['ixon_raw_hist' num2str(kk)]);end
    end

    % Plot raw count total
    raw_opts=struct;
    raw_opts.FitLinear=0;
    hF_ixon_rawtotal=ixon_showRawCountTotal(ixondata,ixon_xVar,raw_opts);

    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_rawtotal,['ixon_raw_counts']);end

end

%% Calculate FFT
ixon_doFFT=0;

fft_opts=struct;
fft_opts.doSmooth=1;
fft_opts.smoothRadius=5;

if ixon_doFFT
    ixondata=ixon_computeFFT(ixondata,fft_opts);
end


%% ANALYSIS : BOX COUNT
ixon_doBoxCount=1;

if ixon_doBoxCount
    ixondata=ixon_boxCount(ixondata);
end

%% PLOTTING : BOX COUNT

ixon_boxPopts=struct;
ixon_boxPopts.NumberExpFit = 0;
ixon_boxPopts.NumberLorentzianFit=0;

ixon_boxPopts.CenterSineFit = 0;       % Fit sine fit to cloud center
ixon_boxPopts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
ixon_boxPopts.CenterLinearFit = 0;     % Linear fit to cloud center

if ixon_doBoxCount  
    % Plot the atom number
    [hF_ixon_numberbox,Ndatabox]=ixon_showBoxNumber(ixondata,ixon_xVar,ixon_boxPopts);      
    yl=get(gca,'YLim');
    set(gca,'YLim',[0 yl(2)]); 
    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_numberbox,'ixon_box_number');end     
    
    % Plot the second moments
    hF_ixon_size=ixon_showBoxMoments(ixondata,ixon_xVar);   
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_size,'ixon_box_size');end     
    
    % Plot the cloud center
    hF_ixon_center=ixon_showBoxCentre(ixondata,ixon_xVar,ixon_boxPopts); 
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_center,'ixon_box_centre');end 
end

%% ANALYSIS : 2D Gaussian
ixon_doGaussFit=0;
% do a very basic PCA to determine angle of the atomic cloud
ixondata=ixon_simple_pca(ixondata);

ixon_gauss_opts=struct;
ixon_gauss_opts.doRescale=1;     % Rescaling the image makes fitting faster
ixon_gauss_opts.doMask=1;        % Apply the image mask
ixon_gauss_opts.Scale=0.5;       % Scale to rescale the image by
ixon_gauss_opts.doRotate=1;      % Allow for gaussian to be rotated (requires PCA)
ixon_gauss_opts.Mask=ixon_mask;  % The image mask

if ixon_doGaussFit  
    ixondata=ixon_gaussFit(ixondata,ixon_gauss_opts);
end
%% PLOTTING : GAUSSIAN
ixon_gauss_opts.NumberExpFit = 0;        % Fit exponential decay to atom number
ixon_gauss_opts.NumberLorentzianFit=0;   % Fit atom number to lorentzian

ixon_gauss_opts.CenterSineFit = 0;       % Fit sine fit to cloud center
ixon_gauss_opts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
ixon_gauss_opts.CenterParabolaFit = 0;
ixon_gauss_opts.CenterLinearFit = 0;     % Linear fit to cloud center

if ixon_doGaussFit
    % Statistics if no variable is changing
    if isequal(ixon_xVar,'ExecutionDate')
        hF_stats=ixon_showGaussStats(ixondata);     
        if ixon_doSave;ixon_saveFigure(ixondata,hF_stats,'ixon_gauss_stats');end
    end
       
    % Counts
    [hF_numbergauss,Ndatagauss]=ixon_showGaussNumber(ixondata,ixon_xVar,ixon_gauss_opts);  
     %ylim([0 max(get(gca,'YLim'))]);
     %ylim([3.5E6 4.5E6]);
     %xlim([0 max(get(gca,'XLim'))]);         
    if ixon_doSave;ixon_saveFigure(ixondata,hF_numbergauss,'ixon_gauss_number');end    
    
    % Size
    hF_size=ixon_showGaussSize(ixondata,ixon_xVar);    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_size,'ixon_gauss_size');end
        
    % Aspect Ratio
    hF_ratio=ixon_showGaussAspectRatio(ixondata,ixon_xVar);    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ratio,'ixon_gauss_ratio');end
    
    % Centre
    hF_Centre=ixon_showGaussCentre(ixondata,ixon_xVar,ixon_gauss_opts);    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_Centre,'ixon_gauss_position');end
        
     % Style of profile --> cut or sum?
    style='cut';
%     style='sum';
    clear hF_X;    
    clear hF_Y;
    hF_X=[];
    hF_Y=[];
    
    hF_Xs=ixon_showGaussProfile(ixondata,'X',style,ixon_xVar);        
    hF_Ys=ixon_showGaussProfile(ixondata,'Y',style,ixon_xVar);  

%   Save the figures (this can be slow)
    if ixon_doSave
        for kk=1:length(hF_Xs)            
            ixon_saveFigure(ixondata,hF_Xs(kk),['ixon_gauss_profile_X' num2str(rNum) '_' num2str(kk)]);
        end
        for kk=1:length(hF_Ys)
            ixon_saveFigure(ixondata,hF_Ys(kk),['ixon_gauss_profile_Y' num2str(rNum) '_' num2str(kk)]);
        end
    end
    
end

%% Animate cloud
ixon_doAnimate = 1;
if ixon_doAnimate == 1
    ixon_animateOpts=struct;
    ixon_animateOpts.StartDelay=2; % Time to hold on first picture
    ixon_animateOpts.MidDelay=.25;     % Time to hold in middle picutres
    ixon_animateOpts.EndDelay=2;     % Time to hold final picture

    % Animate in ascending or descending order?
    % animateOpts.Order='descend';    % Asceneding or descending
    ixon_animateOpts.Order='ascend';
    
    % Color limit for image
    ixon_animateOpts.CLim=[0 500];   % Color limits
%      ixon_animateOpts.CLim='auto';   % Automatically choose CLIM?

    ixon_animate(ixondata,ixon_xVar,ixon_animateOpts);
end

%% Animate cloud FFT
ixon_doAnimateFFT = 1;
if ixon_doAnimateFFT == 1 && ixon_doFFT
    ixon_animateOptsFFT=struct;
    ixon_animateOptsFFT.StartDelay=2; % Time to hold on first picture
    ixon_animateOptsFFT.MidDelay=.25;     % Time to hold in middle picutres
    ixon_animateOptsFFT.EndDelay=2;     % Time to hold final picture

    % Animate in ascending or descending order?
    % animateOpts.Order='descend';    % Asceneding or descending
    ixon_animateOptsFFT.Order='ascend';
    
    % Color limit for image
%     ixon_animateOptsFFT.CLim=[0 1E5];   % Color limits
     ixon_animateOptsFFT.CLim='auto';   % Automatically choose CLIM?

    ixon_animateFFT(ixondata,ixon_xVar,ixon_animateOptsFFT);
end