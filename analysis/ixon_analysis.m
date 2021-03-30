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

global imgdir
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

varType='param';    % The variable is from the sequencer
varType='acq';      % The variable is from the camera

xVar='ExposureTime';
unit='s';


%% Select image directory
% Choose the directory where the images to analyze are stored
disp([datestr(now,13) ' Choose an image analysis folder...']);
dialog_title='Choose the root dire ctory of the images';
imgdir=uigetdir(getImageDir(datevec(now)),dialog_title);
if isequal(imgdir,0)
    disp('Canceling.');
    return 
end

%% Load the data
clear atomdata
disp(['Loading data from ' imgdir]);
files=dir([imgdir filesep '*.mat']);
files={files.name};


for kk=1:length(files)
    str=fullfile(imgdir,files{kk});
    [a,b,c]=fileparts(str);      
    disp(['     (' num2str(kk) ')' files{kk}]);    
    data=load(str);     
    data=data.data;  

    % Display image properties
    try
        disp(['     Image Name     : ' data.Name]);
        disp(['     Execution Time : ' datestr(data.Date)]);
%         disp(['     ' xVar ' : ' num2str(data.Params.(xVar))]);
        disp(' ');
    end    
    
%     if isequal(xVar,'ExecutionDate')
%         data.Params.(xVar)=datenum(data.Params.(xVar))*24*60*60;
%     end
    atomdata(kk)=data;    
end
disp(' ');

if isequal(xVar,'ExecutionDate')
   p=[atomdata.Params] ;
   tmin=min([p.ExecutionDate]);
   for kk=1:length(atomdata)
      atomdata(kk).Params.ExecutionDate= ...
          atomdata(kk).Params.ExecutionDate-tmin;
   end     
end

%% Sort the data
% Sort the data by your given parameter
clear x
disp(['Sorting atomdata by the given ''' xVar '''']);

% Get the object that contains the variable
switch varType
    case 'param'
        varList=[atomdata.Params];
    case 'acq'
        varList=[atomdata.AcquisitionInformation];        
    otherwise
        error('uhh you chose the wrong thing to plot');
end
       
% Make sure that all atomdata have the parameter and record it
allGood=1;
for kk=1:length(varList)
    if isfield(varList(kk),xVar)
        x(kk)=varList.(xVar);
    else
        allGood=0;
        warning(['atomdata(' num2str(kk) ') has no ''' xVar '''']);
    end
end

% Sort if all the data has that parameter
if allGood
    [~, inds]=sort(x);
    atomdata=atomdata(inds);
end

% Get the object that contains the variable
switch varType
    case 'param'
        varList=[atomdata.Params];
    case 'acq'
        varList=[atomdata.AcquisitionInformation];        
    otherwise
        error('uhh you chose the wrong thing to plot');
end
%% Uhh okay

atomdata=computeRawCounts(atomdata);

hist_opts=struct;
hist_opts.GlobalLimits=0;
hist_opts.BinWidth=1;
hist_opts.ImageNumber=1;
hist_opts.YScale='Log';
% hist_opts.YScale='Linear';

for kk=1:size(atomdata(1).RawImages,3)
    disp(kk)
    hist_opts.ImageNumber=kk;
    hF_rawhist=showRawCountHistogram(atomdata,xVar,hist_opts);
    saveFigure(atomdata,hF_rawhist,['raw_hist' num2str(kk)]);
end


raw_opts=struct;
raw_opts.FitLinear=1;
hF_rawtotal=showRawCountTotal(atomdata,xVar,raw_opts);
saveFigure(atomdata,hF_rawtotal,['raw_counts']);
