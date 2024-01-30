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

if ~exist('ixon_auto_dir')
   ixon_auto_dir = 1; 
end

if ixon_auto_dir

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
else
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
ixon_xVar           = 'conductivity_mod_time'; % Variable Name
% ixon_xVar           = 'z_repop_ramptime'; % Variable Name
ixon_overrideUnit   = 'V';    % If ixon_autoUnit=0, use this
ixon_doSave         = 1;    % Save Analysis?
ixon_Magnification = 83;        % Magnification of imaging system
ixon_PixelSize = 16;            % Pixel size in um
%% Analysis Options
% Fitting options
ixon_doBoxCount             = 1;
ixon_doGaussFit             = 0;

% Analysis to run
ixon_doStandardAnalysis     = 1;
ixon_doPlotProfiles         = 0;
ixon_doAnimate              = 0;    % Animate in position domain
ixon_doAnalyzeRaw           = 0;    % Raw Image Analysis
ixon_doAnalyzeFourier       = 0;    % Fourier Domain Analysis
ixon_doAnalyzeStripes2D     = 1;    % Stripe Analysis :  for field stability in titled plane selection

% QGM Single Plane Analysis
ixon_doQGM                  = 1;


%% Image Processing Options

% What do you do to the raw data?
maskname=fullfile('ixon_mask.mat');
ixon_mask=load(maskname);
ixon_mask=ixon_mask.BW;

img_opt = struct;
img_opt.doSubtractBias      = 1;        % Subtract 200 count electronic offset
img_opt.doSubtractBG        = 1;
img_opt.doScale             = 1;        % Scale up image? (good for single-site)
img_opt.ScaleFactor         = 2;        % Amount to scale up by (x2 is good)
img_opt.doRotate            = 1;        % Rotate image? (useful to align along lattices)
img_opt.Theta               = 60.2077;  % Rotation amount (deg.)
img_opt.Theta               = 59.64;  % Rotation amount (deg.)
% img_opt.Theta               = 59.81;  % Rotation amount (deg.)

img_opt.doMask              = 0;        % Mask the data? (not used)
img_opt.Mask                = ixon_mask;% Mask File 512x512
img_opt.doGaussFilter       = 0;        % Filter the image? (bad for single-site)
img_opt.GaussFilterRadius   = 1;        % Filter radius
img_opt.doPSF               = 1;        % Deconvolve with PSF
img_opt.PSF                 = [1.3163 51 21]; % PSF parameters [sigma, N, Niter]
img_opt.doFFT               = 1;        % Compute FFT?
img_opt.doMaskIR            = 1;        % Mask long distance in FFT (useful)
img_opt.IRMaskRadius        = 0.01;     % Mask radius in 1/px
img_opt.doFFTFilter         = 0;        % Filter FFT?
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
[ixondata.Magnification] = deal(ixon_Magnification);
[ixondata.PixelSize] = deal(ixon_PixelSize);

%% Process Images
ixondata = ixonProcessImages(ixondata,img_opt);

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


%% Profiles
profile_opts = struct;
profile_opts.Style = 'cut'; 'sum';  % Cut or sum?
% profile_opts.Style = 'sum';  % Cut or sum?

profile_opts.FigLabel = FigLabel;

clear hF_X;clear hF_Y;
hF_X=[];hF_Y=[];

if ixon_doPlotProfiles

    for rNum=1:size(ixondata(1).ROI,1)
        profile_opts.ROINum = rNum;

        hF_Xs_rNum=ixon_showProfile(ixondata,'X',ixon_xVar,profile_opts);

        if ixon_doSave
            for kk=1:length(hF_Xs_rNum) 
                figure(hF_Xs_rNum(kk));
                ixon_saveFigure2(hF_Xs_rNum(kk),['ixon_R' num2str(rNum) '_X' num2str(kk)],saveOpts);
                pause(0.1);
            end 
        end

        hF_Ys_rNum=ixon_showProfile(ixondata,'Y',ixon_xVar,profile_opts);          
    %   Save the figures (this can be slow)
        if ixon_doSave        
            for kk=1:length(hF_Ys_rNum)
                figure(hF_Ys_rNum(kk));
                ixon_saveFigure2(hF_Ys_rNum(kk),['ixon_R' num2str(rNum) '_Y' num2str(kk)],saveOpts);
                pause(0.1);
            end
        end
        hF_X=[hF_X; hF_Xs_rNum];
        hF_Y=[hF_Y; hF_Ys_rNum];
    end   
end

%% Animate cloud 
if ixon_doAnimate == 1 && ixon_doSave
    ixon_animateOpts=struct;
    
    ixon_animateOpts.xUnit=ixon_unit;
    ixon_animateOpts.StartDelay=2; % Time to hold on first picture
    ixon_animateOpts.MidDelay=.8;     % Time to hold in middle picutres
    ixon_animateOpts.EndDelay=2;     % Time to hold final picture

    % Animate in ascending or descending order?
    % animateOpts.Order='descend';    % Asceneding or descending
    ixon_animateOpts.Order='ascend';
    
    % Color limit for image
ixon_animateOpts.Source = 'ZNoFilter';

     ixon_animateOpts.CLim='auto';   % Automatically choose CLIM?
    ixon_animateOpts.CLim=[0 100];   % Color limits
%     ixon_animateOpts.CLim=[0 1800];   % Color limits


    ixon_animate(ixondata,ixon_xVar,ixon_animateOpts);
end

%% Standard Cloud Analysis
if ixon_doStandardAnalysis;ixon_StandardAnalysis;end
%% Raw Image Analysis
if ixon_doAnalyzeRaw;ixon_AnalyzeRawImages;end
%% Fourier Analysis
if ixon_doAnalyzeFourier;ixon_AnalyzeFourier;end
%% Quantum Gas Micrscopy
if ixon_doQGM; ixon_QGM;end
%% Stripe Analysis
if ixon_doAnalyzeStripes2D;ixon_stripe_2d;end