% ixonAnalysis.m

% This is an imaging analysis script. 
%   Date
disp(repmat('-',1,60));    
disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));    
disp(repmat('-',1,60));    

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))
    

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
global doRotate

lambda=770E-9;

% Choose display rotation
doRotate=0;

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



xVar='ExecutionDate';
unit='s';


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
        disp(['     ' xVar ' : ' num2str(data.Params.(xVar))]);
        disp(' ');
    end    
    
    if isequal(xVar,'ExecutionDate')
        data.Params.(xVar)=datenum(data.Params.(xVar))*24*60*60;
    end
    
    ixondata(kk)=data;    
end
disp(' ');

if isequal(xVar,'ExecutionDate')
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
disp(['Sorting ixondata by the given ''' xVar '''']);

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
    if isfield(varList(kk),xVar)
        x(kk)=varList(kk).(xVar);
    else
        allGood=0;
        warning(['ixondata(' num2str(kk) ') has no ''' xVar '''']);
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
    
    Z=imgs(:,:,2)-imgs(:,:,1);
    ixondata(kk).Z=Z;
end


%% Basic Raw Image Analysis

doRawImageAnalysis=1;
if doRawImageAnalysis   

% Do basic analysis on raw counts
ixondata=ixon_computeRawCounts(ixondata);

% Plot histogram of raw counts
hist_opts=struct;
hist_opts.Outliers=[10 50]; % Histogram wont plot outliers of this many low/high
hist_opts.GlobalLimits=1;   % Maintain historgram x limits
hist_opts.BinWidth=1;       % Histogram bin width
hist_opts.ImageNumber=1;    % Which image to histogram (overwritten)
hist_opts.YScale='Log';     % Histogram y scale
% hist_opts.YScale='Linear';

for kk=1:size(ixondata(1).RawImages,3)
    hist_opts.ImageNumber=kk;
    hF_ixon_rawhist=ixon_showRawCountHistogram(ixondata,xVar,hist_opts);
    if ixon_doSave
        ixon_saveFigure(ixondata,hF_ixon_rawhist,['ixon_raw_hist' num2str(kk)]);
    end
end

% Plot raw count total
raw_opts=struct;
raw_opts.FitLinear=0;
hF_ixon_rawtotal=ixon_showRawCountTotal(ixondata,xVar,raw_opts);

if ixon_doSave
    ixon_saveFigure(ixondata,hF_ixon_rawtotal,['ixon_raw_counts']);
end

end

%% ANALYSIS : BOX COUNT
doBoxCount=1;

if doBoxCount
    ixondata=ixon_boxCount(ixondata);
end


%% PLOTTING : BOX COUNT

ixon_boxPopts = struct;
ixon_boxPopts.NumberExpFit = 0;                  % Fit exponential decay to atom number

ixon_boxPopts.CenterSineFit = 0;       % Fit sine fit to cloud center
ixon_boxPopts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
ixon_boxPopts.CenterLinearFit = 0;     % Linear fit to cloud center

if doBoxCount  
    % Plot the atom number
    [hF_ixon_numberbox,Ndatabox]=ixon_showBoxNumber(ixondata,xVar,ixon_boxPopts);      
    yl=get(gca,'YLim');
    set(gca,'YLim',[0 yl(2)]); 
    if ixon_doSave
        ixon_saveFigure(ixondata,hF_ixon_numberbox,'ixon_box_number'); 
    end     
    
    % Plot the second moments
    hF_ixon_size=ixon_showBoxMoments(ixondata,xVar);   
    if ixon_doSave
        ixon_saveFigure(ixondata,hF_ixon_size,'ixon_box_size'); 
    end     
    
    % Plot the cloud center
    hF_ixon_center=ixon_showBoxCentre(ixondata,xVar,ixon_boxPopts); 
    if ixon_doSave
        ixon_saveFigure(ixondata,hF_ixon_center,'ixon_box_centre'); 
    end 
end

%% Animate cloud
ixon_doAnimate = 1;
if ixon_doAnimate == 1
ixon_animateOpts=struct;
ixon_animateOpts.StartDelay=2; % Time to hold on first picture
ixon_animateOpts.MidDelay=.25;     % Time to hold in middle picutres
ixon_animateOpts.EndDelay=2;     % Time to hold final picture


% animateOpts.Order='descend';    % Asceneding or descending
ixon_animateOpts.Order='ascend';
ixon_animateOpts.CLim=[0 20000];   % Color limits

ixon_animate(ixondata,xVar,ixon_animateOpts);
end
