function ixondata = ixon_main_focus
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
end

imgdir = newdir;
ixon_imgdir = imgdir;
saveDir = [imgdir filesep 'figures'];
if ~exist(saveDir,'dir'); mkdir(saveDir);end   
saveOpts.saveDir=saveDir;
strs=strsplit(imgdir,filesep);
FigLabel=[strs{end-1} filesep strs{end}];


%% Analysis Variable
% This section of code chooses the variable to plot against for aggregate
% plots.  The chosen variable MUST match a variable provided in the params
% field of the .mat file. The unit has no tangibile affect and only affects
% display properties.

% Choose what kind of variable to plot against (sequencer/camera)
varType             = 'param'; % always select 'param' for now 
ixon_autoXVar       = 0;      % Auto detect changing variable?
ixon_autoUnit       = 1;      % Auto detect unit for variable?
ixon_xVar           = 'ExecutionDate'; % Variable Name
ixon_overrideUnit   = 'V';    % If ixon_autoUnit=0, use this
ixon_doSave         = 1;    % Save Analysis?
ixon_Magnification  = 83;        % Magnification of imaging system
ixon_PixelSize      = 16;            % Pixel size in um


% Ignore these variables when choosing auto var
autoVar_Ignore = {'f_offset','piezo_offset'};
% autoVar_Ignore = {};

%% Analysis Options

ixon_doBoxCount                     = 1;

%% Image Processing Options

% What do you do to the raw data?
maskname=fullfile('ixon_mask.mat');
ixon_mask=load(maskname);
ixon_mask=ixon_mask.BW;

img_opt = struct;
img_opt.doSubtractBias      = 1;        % Subtract 200 count electronic offset
img_opt.doSubtractBG        = 1;
img_opt.doScale             = 0;        % Scale up image? (good for single-site)
img_opt.ScaleFactor         = 2;        % Amount to scale up by (x2 is good)
img_opt.doRotate            = 0;        % Rotate image? (useful to align along lattices)
img_opt.Theta               = 59.81;  % Rotation amount (deg.)
img_opt.DetectNoise         = 1;
img_opt.doMask              = 0;        % Mask the data? (not used)
img_opt.Mask                = ixon_mask;% Mask File 512x512
img_opt.doGaussFilter       = 0;        % Filter the image? (bad for single-site)
img_opt.GaussFilterRadius   = 1;        % Filter radius
img_opt.doPSF               = 0;    % Deconvolve with PSF
img_opt.PSF                 = [1.3163 11 51]; % PSF parameters [sigma, N, Niter]
img_opt.doFFT               = 0;        % Compute FFT?
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
    
    for kk=1:length(autoVar_Ignore)
        thisVar = autoVar_Ignore{kk};
        index = find(ismember(xVars, thisVar));
        xVars(index)=[];
    end
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
ixondata = ixon_ProcessImages(ixondata,img_opt);

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


%% Focus Analysis

 focus = ixon_MultiShotFocus(ixondata);

[hF]= ixon_showMultiShotFocusSummary(focus,ixon_xVar);




end

