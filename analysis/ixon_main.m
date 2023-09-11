% ixonAnalysis.m
%
% Author : CF Fujiwara
%
% This script is the primary analysis file for the iXon MATLAB code. It
% calls and plots all other analyses.

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));    

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

%% Select image directory
% Choose the directory where the images to analyze are stored
disp([datestr(now,13) ' Choose an image analysis folder...']);
dialog_title='Choose the root directory of the images';
default_analysis_dir = ixon_getDayDir;
saveOpts = struct;
saveOpts.Quality = 'auto';

newdir=uigetdir(default_analysis_dir,dialog_title);
if isequal(newdir,0)
    disp('Canceling.');    
    return; 
else
    imgdir = newdir;
    ixon_imgdir = imgdir;
    saveDir = [imgdir filesep 'figures'];
    if ~exist(saveDir,'dir'); mkdir(saveDir);end   
    saveOpts.saveDir=saveDir;
    strs=strsplit(imgdir,filesep);
    FigLabel=[strs{end-1} filesep strs{end}];
end

%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% Choose what kind of variable to plot against (sequencer/camera)
varType             = 'param'; % always select 'param' for now 
ixon_autoXVar       = 1;      % Auto detect changing variable?
ixon_autoUnit       = 1;      % Auto detect unit for variable?
ixon_xVar           = 'qgm_raman_2photon_detuning'; % Variable Name
ixon_overrideUnit   = 'V';    % If ixon_autoUnit=0, use this
ixon_doSave         = 1;    % Save Analysis?

%% Analysis Options
% Select what kinds of analyses you'd like to perform
doRawImageHistogram=0;
ixon_doBoxCount=1;
ixon_doGaussFit=0;

% Fast Fourier Transform Analysis
% Use if you are looking for astigmatism in the image
ixon_doFFT=0;
ixon_fft_doBoxCount=0;

% Stripe Analysis
% This is used to analyze the field during the plane selection
do_2dStripeAnalysis=0;
doStripeAnalysis=0;

ixon_doAnimate = 1;
%% Image Processing Options

% What do you do to the raw data?
maskname=fullfile('ixon_mask.mat');
ixon_mask=load(maskname);
ixon_mask=ixon_mask.BW;

img_opt = struct;
img_opt.doSubtractBias      = 1;        % Subtract 200 count electronic offset
img_opt.doScale             = 1;        % Scale up image? (good for single-site)
img_opt.ScaleFactor         = 2;        % Amount to scale up by (x2 is good)
img_opt.doRotate            = 1;        % Rotate image? (useful to align along lattices)
img_opt.Theta               = 60.2077;  % Rotation amount (deg.)
img_opt.doMask              = 0;        % Mask the data? (not used)
img_opt.Mask                = ixon_mask;% Mask File 512x512
img_opt.doGaussFilter       = 0;        % Filter the image? (bad for single-site)
img_opt.GaussFilterRadius   = 1;        % Filter radius
img_opt.doPSF               = 0;        % Deconolve with PSF
img_opt.PSF                 = [1.3163 51 12]; % PSF parameters [sigma, N, Niter]
img_opt.doFFT               = 1;        % Compute FFT?
img_opt.doMaskIR            = 1;        % Mask long distance in FFT (useful)
img_opt.IRMaskRadius        = 0.01;     % Mask radius in 1/px
img_opt.doFFTFilter         = 1;        % Filter FFT?
img_opt.FFTFilterRadius     = 1;        % FFT Filter radius (1/px)

%% Analysis ROI
% Analysis ROI is an Nx4 matrix of [X1 X2 Y1 Y2] which specifies a region
% to analyze. Each new row in the matrix indicates a separate ROI to
% perform analysis on.

% Full ROI 
ixonROI = [1 512 1 512]; 

%% Load the data
clear ixondata
disp(['Loading data from ' ixon_imgdir]);
files=dir([ixon_imgdir filesep '*.mat']);
files={files.name};

for kk=1:length(files)
    str=fullfile(ixon_imgdir,files{kk});
    [a,b,c]=fileparts(str);      
    disp(['     (' num2str(kk) ')' files{kk}]);    
    data=load(str);data=data.data;  
    try
        disp(['     Image Name     : ' data.Name]);
        disp(['     Execution Time : ' datestr(data.Date)]);
        if ~ixon_autoXVar
            disp(['     ' ixon_xVar ' : ' num2str(data.Params.(ixon_xVar))]);
        end
        disp(' ');
    end    
    data.Params.ExecutionDate = datenum(data.Params.ExecutionDate);
    data.Params.ExecutionDateStr = datestr(data.Params.ExecutionDate); 
    ixondata(kk)=data;    
end
disp(' ');
%% Match Parameters and Flags
% This makes sure that all data as has all the flags and params in each. 
% This is usually not necessary.
ixondata = ixon_matchParamsFlags(ixondata);

%% X Variable and Units
% If auto unit and variable are chosen, search through the parameters and
% data to find which variable(s) are being changed.

% Also, sort the data by that chosen variable.

if ixon_autoXVar
    xVars = ixon_findXVars(ixondata);
    disp([' Found ' num2str(length(xVars)) ...
        ' valid variables that are changing to plot against.']);
    disp(xVars);    
    % Select the first variable that is changed
    ind = 1;    
    ixon_xVar = xVars{ind};    
    disp([' Setting ' ixon_xVar ' to be the x-variable']);    
    for kk=1:length(ixondata)
        disp([' (' num2str(kk) ') (' num2str(ixondata(kk).Params.(ixon_xVar)) ') ' ...
            ixondata(kk).Name]); 
    end
    disp(' ');
end

if ixon_autoUnit && isfield(ixondata(1),'Units')  && isequal(varType,'param')
    ixon_unit=ixondata(1).Units.(ixon_xVar);
else
    ixon_unit=ixon_overrideUnit;
end

if isequal(ixon_xVar,'ExecutionDate')
   ixon_unit='days'; 
end

% Sort the data by your given parameter
disp(['Sorting atomdata by the given ''' ixon_xVar '''']);
x=zeros(length(ixondata),1);
for kk=1:length(ixondata)
    if isfield(ixondata(kk).Params,ixon_xVar)
        x(kk)=ixondata(kk).Params.(ixon_xVar);
    else
        warning(['atomdata(' num2str(kk) ') has no ''' ixon_xVar '''']);
    end
end

% Sort it
[~, inds]=sort(x);
ixondata=ixondata(inds);

% Get the object that contains the variable
switch varType
    case 'param'
        varList=[ixondata.Params];
    case 'acq'
        varList=[ixondata.AcquisitionInformation];        
    otherwise
        error('uhh you chose the wrong thing to plot');
end
%% Save Params
if ixon_doSave
    Params =[ixondata.Params];
    filename=fullfile(ixon_imgdir,'figures','Params.mat');
    save(filename,'Params');
end
%% Distribute ROI
[ixondata.ROI]=deal(ixonROI);
%% Process Images
ixondata = ixonProcessImages(ixondata,img_opt);
%% Basic Raw Image Analysis

if doRawImageHistogram   
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
%% ANALYSIS : BOX COUNT

if ixon_doBoxCount
    ixondata=ixon_boxCount(ixondata);    
    ixon_boxdata = ixon_getBoxData(ixondata,ixon_xVar);
    ixon_boxdata.ProcessOptions = img_opt;
    if ixon_doSave
        P = [ixon_boxdata.Params];
        BoxData_SourceFiles = {P.ExecutionDateStr};
        BoxData_xVarName = ixon_xVar;
        BoxData_xVarUnit = ixon_unit;
        BoxData_xVar = [P.(ixon_xVar)];
        BoxData_N = [ixon_boxdata.N];
        BoxData_Xc = [ixon_boxdata.Xc];
        BoxData_Yc = [ixon_boxdata.Yc];
        BoxData_Xs = [ixon_boxdata.Xs];
        BoxData_Ys = [ixon_boxdata.Ys];

        filename=fullfile(ixon_imgdir,'figures','ixon_boxdata_python.mat');
        save(filename,'BoxData_SourceFiles','BoxData_xVarName','BoxData_xVarUnit',...
            'BoxData_xVar','BoxData_N','BoxData_Xc','BoxData_Yc',...
            'BoxData_Xs','BoxData_Ys');
        filename=fullfile(ixon_imgdir,'figures','ixon_boxdata.mat');
        save(filename,'ixon_boxdata');
    end
end

%% ANALYSIS : 2D Gaussian
% do a very basic PCA to determine angle of the atomic cloud
if ixon_doGaussFit  
    ixon_gauss_opts=struct;
    ixon_gauss_opts.doRescale=1;     % Rescaling the image makes fitting faster
    ixon_gauss_opts.doMask=1;        % Apply the image mask
    ixon_gauss_opts.Scale=0.25;       % Scale to rescale the image by
    ixon_gauss_opts.doRotate=1;      % Allow for gaussian to be rotated (requires PCA)
    ixon_gauss_opts.Mask=ixon_mask;  % The image mask
    ixon_gauss_opts.doBackground = 1; % Enable a background to the fit
    
    ixondata=ixon_gaussFit(ixondata,ixon_gauss_opts);
    ixon_gaussdata = ixon_getGaussData(ixondata,ixon_xVar);
    ixon_gaussdata.ProcessOptions = img_opt;    
    
    if ixon_doSave
        P = [ixon_gaussdata.Params];
        GaussData_SourceFiles = {P.ExecutionDateStr};
        GaussData_xVarName = ixon_xVar;
        GaussData_xVarUnit = ixon_unit;
        GaussData_xVar = [P.(ixon_xVar)];
        GaussData_N = [ixon_gaussdata.N];
        GaussData_Xc = [ixon_gaussdata.Xc];
        GaussData_Yc = [ixon_gaussdata.Yc];
        GaussData_Xs = [ixon_gaussdata.Xs];
        GaussData_Ys = [ixon_gaussdata.Ys];
        filename=fullfile(ixon_imgdir,'figures','ixon_gaussdata_python.mat');
        save(filename,'GaussData_SourceFiles','GaussData_xVarName','GaussData_xVarUnit',...
            'GaussData_xVar','GaussData_N','GaussData_Xc','GaussData_Yc',...
            'GaussData_Xs','GaussData_Ys');
        filename=fullfile(ixon_imgdir,'figures','ixon_gaussdata.mat');
        save(filename,'ixon_gaussdata');
    end    
end

%% Calculate FFT
% 
% fft_opts=struct;
% fft_opts.doSmooth=1;
% fft_opts.smoothRadius=1;
% fft_opts.fft_N=2^11; % Can go higher for smoother data
% 
% fft_opts.maskIR=0;
% fft_opts.maskUV=0;
% fft_opts.LMax=50;
% fft_opts.LMin=5;
% 
% 
% if ixon_doFFT
%     ixondata=ixon_computeFFT(ixondata,fft_opts);
% end
% 
% % Apply makss to FFT Data
% ixon_mask_IR=1;
% ixon_mask_UV=0;
% 
% 
% if fft_opts.maskIR
%     ixondata=ixon_fft_maskIR(ixondata,fft_opts.LMax);    
% end
% 
% if fft_opts.maskUV
%     ixondata=ixon_fft_maskUV(ixondata,fft_opts.LMin);    
% end

%% ANALYSIS : FFT BOX COUNT

% if ixon_fft_doBoxCount && ixon_doFFT    
%     fft_boxOpts=struct;
%     fft_boxOpts.maskIR=fft_opts.maskIR;
%     fft_boxOpts.LMax=fft_opts.LMax;
%     fft_boxOpts.maskUV=fft_opts.maskUV;
%     fft_boxOpts.LMin=fft_opts.LMin;
%         ixondata=ixon_fft_boxCount(ixondata,fft_boxOpts);
% end

%% PLOTTING : FFT BOX COUNT
% ixon_fft_boxPopts=struct;
% ixon_fft_boxPopts.xUnit=ixon_unit;
% if ixon_fft_doBoxCount  && ixon_doFFT 
%     % Plot the second moments
%     hF_ixon_size=ixon_fft_showBoxMoments(ixondata,ixon_xVar,ixon_fft_boxPopts);   
%     if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_size,'ixon_fft_box_size');end          
% end

%% PLOTTING : BOX COUNT

ixon_boxPopts=struct;
ixon_boxPopts.xUnit=ixon_unit;

% ixon_boxPopts.NumberScale='Log';
ixon_boxPopts.NumberScale='Linear';

ixon_boxPopts.NumberExpFit = 0;
ixon_boxPopts.NumberExp2SumFit = 0;

ixon_boxPopts.NumberLorentzianFit=0;

ixon_boxPopts.CenterSineFit = 0;       % Fit sine fit to cloud center
ixon_boxPopts.CenterDecaySineFit = 1;  % Fit decaying sine to cloud center
ixon_boxPopts.CenterGrowSineFit = 0;  % Fit decaying sine to cloud center
ixon_boxPopts.CenterLinearFit = 0;     % Linear fit to cloud center

if ixon_doBoxCount  
    % Plot the atom number
    [hF_ixon_numberbox,Ndatabox]=ixon_showBoxNumber(ixondata,ixon_xVar,ixon_boxPopts);      
    yl=get(gca,'YLim');
    set(gca,'YLim',[0 yl(2)]);
%     set(gca,'YLim',[2.0e8 2.5e8]);
%     set(gca,'XScale','log');
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_numberbox,'ixon_box_number');end     
    
    % Plot the second moments
    hF_ixon_size=ixon_showBoxMoments(ixondata,ixon_xVar,ixon_boxPopts);   
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_size,'ixon_box_size');end     
    
    % Plot the cloud center
    [hF_ixon_center,Xc,Yc]=ixon_showBoxCentre(ixondata,ixon_xVar,ixon_boxPopts); 
    if ixon_doSave;ixon_saveFigure(ixondata,hF_ixon_center,'ixon_box_centre');end 

    
end



%% PLOTTING : GAUSSIAN
ixon_gauss_opts.xUnit=ixon_unit;
ixon_gauss_opts.NumberExpFit = 0;        % Fit exponential decay to atom number
ixon_gauss_opts.NumberLorentzianFit=0;   % Fit atom number to lorentzian
ixon_gauss_opts.NumberScale = 'linear'; 
% ixon_gauss_opts.NumberScale = 'log'; 

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
%  xlim([0 1.4]);
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
    
%     hF_stats=ixon_showGaussStats(ixondata,ixon_gauss_opts);     
%     if ixon_doSave;ixon_saveFigure(ixondata,hF_stats,'ixon_stats');end

        
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
    
end


%% Animate cloud 
if ixon_doAnimate == 1 && ixon_doSave
    ixon_animateOpts=struct;
    
    ixon_animateOpts.xUnit=ixon_unit;
    ixon_animateOpts.StartDelay=2; % Time to hold on first picture
    ixon_animateOpts.MidDelay=.1;     % Time to hold in middle picutres
    ixon_animateOpts.EndDelay=2;     % Time to hold final picture

    % Animate in ascending or descending order?
    % animateOpts.Order='descend';    % Asceneding or descending
    ixon_animateOpts.Order='ascend';
    
    % Color limit for image
%     ixon_animateOpts.CLim=[50 200];   % Color limits
%          ixon_animateOpts.CLim=[0 300];   % Color limits

     ixon_animateOpts.CLim='auto';   % Automatically choose CLIM?
%        ixon_animateOpts.CLim=[0 1000];   % Color limits

    ixon_animate(ixondata,ixon_xVar,ixon_animateOpts);
end

%% Animate cloud FFT
% ixon_doAnimateFFT = 1;
% 
% ixon_animateOptsFFT=struct;   
% 
% 
% % Variable to animate versus
% ixon_animateOptsFFT.xUnit=ixon_unit;
% 
% % Animation Timings
% ixon_animateOptsFFT.StartDelay=2; % Time to hold on first picture
% ixon_animateOptsFFT.MidDelay=.25;     % Time to hold in middle picutres
% ixon_animateOptsFFT.EndDelay=2;     % Time to hold final picture
% 
% % Animate in ascending or descending order?
% % animateOpts.Order='descend';    
% ixon_animateOptsFFT.Order='ascend';
% 
% % Color limit for image
% ixon_animateOptsFFT.CLim=[0 .5];   % Color limits 
% ixon_animateOptsFFT.CLim=[0 3];   % Color limits 
% 
% ixon_animateOptsFFT.CLim='auto';   % Automatically choose CLIM?
% 
% % FFT UV Cutoff
% % Reduce the animation view to within a frequency of 1/L
% ixon_animateOptsFFT.mask_UV=1;
% ixon_animateOptsFFT.LMin=20;
% 
% % FFT IR Cutoff
% % Apply mask to interior regions to mask 
% ixon_animateOptsFFT.mask_IR=1;
% ixon_animateOptsFFT.LMax=200;
% 
% if ixon_doAnimateFFT == 1 && ixon_doFFT && ixon_doSave
%     ixon_animateFFT(ixondata,ixon_xVar,ixon_animateOptsFFT);
% end

%% Stripe Analysis
if do_2dStripeAnalysis
   ixon_stripe_2d; 
end

if doStripeAnalysis
   ixon_stripe_1d; 
end

