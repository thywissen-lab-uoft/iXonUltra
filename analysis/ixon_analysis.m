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

% The variable to plot against
ixon_xVar='Objective_Piezo_Z';

% Should the analysis attempt to automatically find the unit?
ixon_autoUnit=1;

% If ixon_autoUnit=0, this will be used.
ixon_overrideUnit='V';

% Flag whether to save the output figures or not (code is faster if not
% saving)
ixon_doSave=1;

% Define the output data
outdata=struct;

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

% if isequal(ixon_xVar,'ExecutionDate')
%    p=[ixondata.Params] ;
%    tmin=min([p.ExecutionDate]);
%    for kk=1:length(ixondata)
%       ixondata(kk).Params.ExecutionDate= ...
%           ixondata(kk).Params.ExecutionDate-tmin;
%    end     
% end


%% Grab the Unit
if ixon_autoUnit && isfield(ixondata(1),'Units')  && isequal(varType,'param')
    ixon_unit=ixondata(1).Units.(ixon_xVar);
else
    ixon_unit=ixon_overrideUnit;
end

if isequal(ixon_xVar,'ExecutionDate')
   ixon_unit='s'; 
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
%% Assign Params to outdata

outdata.Params=ixondata.Params;


%% Analysis ROI
% Analysis ROI is an Nx4 matrix of [X1 X2 Y1 Y2] which specifies a region
% to analyze. Each new row in the matrix indicates a separate ROI to
% perform analysis on.
%
% While in principle different images can have different analysis ROIs,
% this is currently disabled because it creates code issues at the moment.

% Full ROI
ixonROI = [210 297 324 188]; 
ixonROI = [220 324 188 297]; 
ixonROI = [1 512 1 512]; 


[ixondata.ROI]=deal(ixonROI);

%% Image Processing

% Electronic bias offset
doSubBias=1;
offset=200; % (Doesnt affect calculations)

% Ixon mask
doApplyMask=0;
maskname=fullfile('ixon_mask.mat');
ixon_mask=load(maskname);
ixon_mask=ixon_mask.BW;

% Gauss Filter
doGaussFilter=1;
filter_radius=.5;  % Gaussian filter radius


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

%% PSF Sharpening

doSharpenPSF=0;
if doSharpenPSF
   ixondata=sharpenPSF(ixondata); 
end

%% Basic Raw Image Analysis

doRawImageAnalysis=0;
if doRawImageAnalysis   

    % Do basic analysis on raw counts
    ixondata=ixon_computeRawCounts(ixondata);

    % Plot histogram of raw counts
    hist_opts=struct;    
    hist_opts.xUnit=ixon_unit;

    % Specify the plot variable and units    
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
    raw_opts.xUnit=ixon_unit;
    
    % Define the variable and units    
    raw_opts.FitLinear=0;
    
    hF_ixon_rawtotal=ixon_showRawCountTotal(ixondata,ixon_xVar,raw_opts);

    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_rawtotal,['ixon_raw_counts']);end

end

%% Calculate FFT
ixon_doFFT=1;

fft_opts=struct;
fft_opts.doSmooth=1;
fft_opts.smoothRadius=1;
fft_opts.fft_N=2^11; % Can go higher for smoother data

fft_opts.maskIR=0;
fft_opts.maskUV=0;
fft_opts.LMax=50;
fft_opts.LMin=5;


if ixon_doFFT
    ixondata=ixon_computeFFT(ixondata,fft_opts);
end

% Apply makss to FFT Data
ixon_mask_IR=1;
ixon_mask_UV=0;


if fft_opts.maskIR
    ixondata=ixon_fft_maskIR(ixondata,fft_opts.LMax);    
end

if fft_opts.maskUV
    ixondata=ixon_fft_maskUV(ixondata,fft_opts.LMin);    
end

%% ANALYSIS : FFT BOX COUNT
ixon_fft_doBoxCount=0;

fft_boxOpts=struct;
fft_boxOpts.maskIR=fft_opts.maskIR;
fft_boxOpts.LMax=fft_opts.LMax;
fft_boxOpts.maskUV=fft_opts.maskUV;
fft_boxOpts.LMin=fft_opts.LMin;

if ixon_fft_doBoxCount && ixon_doFFT
    ixondata=ixon_fft_boxCount(ixondata,fft_boxOpts);

end

%% PLOTTING : FFT BOX COUNT

ixon_fft_boxPopts=struct;
ixon_fft_boxPopts.xUnit=ixon_unit;

if ixon_fft_doBoxCount  && ixon_doFFT
 
    % Plot the second moments
    hF_ixon_size=ixon_fft_showBoxMoments(ixondata,ixon_xVar,ixon_fft_boxPopts);   
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_size,'ixon_fft_box_size');end     
     
end

%% ANALYSIS : BOX COUNT
ixon_doBoxCount=1;

if ixon_doBoxCount
    ixondata=ixon_boxCount(ixondata);
end

%% PLOTTING : BOX COUNT

ixon_boxPopts=struct;
ixon_boxPopts.xUnit=ixon_unit;

% ixon_boxPopts.NumberScale='Log';
ixon_boxPopts.NumberScale='Linear';

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
%     set(gca,'YLim',[2.0e8 2.5e8]);
    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_numberbox,'ixon_box_number');end     
    
    % Plot the second moments
    hF_ixon_size=ixon_showBoxMoments(ixondata,ixon_xVar,ixon_boxPopts);   
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_size,'ixon_box_size');end     
    
    % Plot the cloud center
    [hF_ixon_center,Xc,Yc]=ixon_showBoxCentre(ixondata,ixon_xVar,ixon_boxPopts); 
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_center,'ixon_box_centre');end 

    outdata.Ndatabox=Ndatabox;
end



%% ANALYSIS : 2D Gaussian
ixon_doGaussFit=1;
% do a very basic PCA to determine angle of the atomic cloud
% ixondata=ixon_simple_pca(ixondata);

ixon_gauss_opts=struct;
ixon_gauss_opts.doRescale=1;     % Rescaling the image makes fitting faster
ixon_gauss_opts.doMask=1;        % Apply the image mask
ixon_gauss_opts.Scale=0.5;       % Scale to rescale the image by
ixon_gauss_opts.doRotate=1;      % Allow for gaussian to be rotated (requires PCA)
ixon_gauss_opts.Mask=ixon_mask;  % The image mask
% ixon_gauss_opts.doBackground = 0; % Enable a background to the fit
if ixon_doGaussFit  
    ixondata=ixon_gaussFit(ixondata,ixon_gauss_opts);
end
%% PLOTTING : GAUSSIAN
ixon_gauss_opts.xUnit=ixon_unit;
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
    hF_size=ixon_showGaussSize(ixondata,ixon_xVar,ixon_gauss_opts);    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_size,'ixon_gauss_size');end
        
    % Aspect Ratio
    hF_ratio=ixon_showGaussAspectRatio(ixondata,ixon_xVar,ixon_gauss_opts);    
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
    
    hF_Xs=ixon_showGaussProfile(ixondata,'X',style,ixon_xVar,ixon_gauss_opts);        
    hF_Ys=ixon_showGaussProfile(ixondata,'Y',style,ixon_xVar,ixon_gauss_opts);  

%   Save the figures (this can be slow)
    if ixon_doSave
        for kk=1:length(hF_Xs)            
            ixon_saveFigure(ixondata,hF_Xs(kk),['ixon_gauss_profile_X'  num2str(kk)]);
        end
        for kk=1:length(hF_Ys)
            ixon_saveFigure(ixondata,hF_Ys(kk),['ixon_gauss_profile_Y' '_' num2str(kk)]);
        end
    end
    
    outdata.Ndatagauss=Ndatagauss;    
end

%% 1D STRIPE ANALYSIS
% This analyzes the stripes for field stability.  This is a 1D analysis,
% where a 1D sine wave with a gaussain envelope is fitted over a sub ROI of
% the data.
%
% This one dimensional analysis requires a variety of inputs to fit.  See
% the "1D Fit Parameters" below.
doStripeAnalysis=0;

stripe_opts=struct;

% X plot unit
strope_opts.xUnit=ixon_unit;

% 1D Fit Parameters
stripe_opts.theta=57;               % Rotation Angle
stripe_opts.rotrange=[220 300];     % Sub region to inspect
stripe_opts.FitType='Sine';         % Fit Type
stripe_opts.LowThreshold=0.2;       % Low ampltude to ignore
stripe_opts.L0=80;                 % Guess wavelength (pixels)
stripe_opts.phi0=pi/2;             % Guess phase (radians)
stripe_opts.B0=0.4;                 % Guess modulation

% Animate
stripe_opts.saveAnimation=1;        % save the animation?
stripe_opts.StartDelay=.25;
stripe_opts.MidDelay=.25;
stripe_opts.EndDelay=.25;

% Field Analysis
field_opts=struct;
field_opts.xUnit=ixon_unit;
field_opts.FieldGradient=210;   % In G/cm
field_opts.LatticeSpacing=532E-9; % in meter
field_opts.FitType='Exp';

if doStripeAnalysis
    [hF_stripe,stripe_data]=analyzeStripes(ixondata,ixon_xVar,stripe_opts);
    
    if ixon_doSave;ixon_saveFigure(ixondata,hF_stripe,'ixon_stripe');end

    field_gradient=210; % field gradient in G/cm
    stripe_data.grabMagnetometer=1;
    stripe_data.Nsmooth=10;    
    
    [hF_field_stripe1,hF_field_stripe2,hF_field_sense]=fieldAnalysis(stripe_data,field_opts);
    
    if ixon_doSave
        ixon_saveFigure(ixondata,hF_field_stripe1,'ixon_field_stripe1');  
        ixon_saveFigure(ixondata,hF_field_stripe2,'ixon_field_stripe2');        

         if stripe_data.grabMagnetometer
            ixon_saveFigure(ixondata,hF_field_sense,'ixon_field_sense');
        end  
    end
    
    outdata.stripe_data=stripe_data;   
end

%% 2D Stripe Analysis
% This analzyes the stripes for take from the ixon camera.  This is a 2D
% fit over the entire cloud.  This fits a 2D gaussian modulated by a sine
% wave at a particular angle.  This is pariticularly useful to fit the
% angular dependence of the data. 
%
% The input fit parameters are specified in the options structure.

do_2dStripeAnalysis=0;

stripe_2d_opts=struct;

stripe_2d_opts.xUnit=ixon_unit;

stripe_2d_opts.ShimFit=0;
stripe_2d_opts.Theta=[-90 90]; % Specify the domain (MUST BE 180 DEGREES)
stripe_2d_opts.saveAnimation=1;        % save the animation?
stripe_2d_opts.StartDelay=1;
stripe_2d_opts.MidDelay=.5;
stripe_2d_opts.EndDelay=1;

if do_2dStripeAnalysis
    [hF_stripe_2d,stripe_data2d]=analyzeStripes2(ixondata,ixon_xVar,stripe_2d_opts);

    if ixon_doSave
        ixon_saveFigure(ixondata,hF_stripe_2d,'ixon_field_stripe_2d');        
    end
    
    outdata.stripe_data2d=stripe_data2d;   
end
%% Animate cloud 
ixon_doAnimate = 1;
if ixon_doAnimate == 1 && ixon_doSave
    ixon_animateOpts=struct;
    
    ixon_animateOpts.xUnit=ixon_unit;
    ixon_animateOpts.StartDelay=2; % Time to hold on first picture
    ixon_animateOpts.MidDelay=1;     % Time to hold in middle picutres
    ixon_animateOpts.EndDelay=2;     % Time to hold final picture

    % Animate in ascending or descending order?
    % animateOpts.Order='descend';    % Asceneding or descending
    ixon_animateOpts.Order='ascend';
    
    % Color limit for image
%     ixon_animateOpts.CLim=[50 200];   % Color limits
%         ixon_animateOpts.CLim=[50 500];   % Color limits

     ixon_animateOpts.CLim='auto';   % Automatically choose CLIM?
%         ixon_animateOpts.CLim=[2000 18000];   % Color limits

    ixon_animate(ixondata,ixon_xVar,ixon_animateOpts);
end

%% Animate cloud FFT
ixon_doAnimateFFT = 1;

ixon_animateOptsFFT=struct;   

% Variable to animate versus
ixon_animateOptsFFT.xUnit=ixon_unit;

% Animation Timings
ixon_animateOptsFFT.StartDelay=2; % Time to hold on first picture
ixon_animateOptsFFT.MidDelay=.25;     % Time to hold in middle picutres
ixon_animateOptsFFT.EndDelay=2;     % Time to hold final picture

% Animate in ascending or descending order?
% animateOpts.Order='descend';    
ixon_animateOptsFFT.Order='ascend';

% Color limit for image
ixon_animateOptsFFT.CLim=[0 .5];   % Color limits 
% ixon_animateOptsFFT.CLim='auto';   % Automatically choose CLIM?

% FFT UV Cutoff
% Reduce the animation view to within a frequency of 1/L
ixon_animateOptsFFT.mask_UV=1;
ixon_animateOptsFFT.LMin=5;

% FFT IR Cutoff
% Apply mask to interior regions to mask 
ixon_animateOptsFFT.mask_IR=1;
ixon_animateOptsFFT.LMax=200;

if ixon_doAnimateFFT == 1 && ixon_doFFT && ixon_doSave
    ixon_animateFFT(ixondata,ixon_xVar,ixon_animateOptsFFT);
end

%% save output data
if ixon_doSave
    filename=fullfile(ixon_imgdir,'figures','outdata.mat');
    save(filename,'outdata');
end
