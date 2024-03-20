function pd_analysis_wrapper(arg1,pd_opts)

if nargin <2
   pd_opts = struct; 
end

if nargin < 1    
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
        ixon_imgdir = newdir;
        saveDir = [ixon_imgdir filesep 'figures'];
        if ~exist(saveDir,'dir'); mkdir(saveDir);end   
        strs=strsplit(ixon_imgdir,filesep);
        FigLabel=[strs{end-1} filesep strs{end}];        
        pd_opts.saveDir = saveDir;
        pd_opts.FigLabel = FigLabel;
    end    
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

%%

params = [ixondata.Params];
names = {ixondata.Name};
dates = {ixondata.Date};
pdsrc = 'X:\LabJackLogs\ODTQPD'; %Windows

for jj=1:length(params)

    % Date string for run
    fdate = datestr(params(jj).ExecutionDateStr,'YYYY-mm-dd_HH-MM-SS');
    
    % Directory of QPD file
    pddir = fullfile(pdsrc,fdate(1:4),[fdate(1:4) '.' fdate(6:7)],[fdate(6:7) '.' fdate(9:10)]);
    
    % QPD filename
    fname = ['ODTQPD_' fdate '.mat'];
    qpdfiles{jj} = fullfile(pddir,fname);
end

%% Load Files
clear pd_data
disp('Loading QPD files');
for nn=1:length(qpdfiles)
    fprintf(['(' num2str(nn) '/' num2str(length(qpdfiles)) ') ']);
    [~,b,~]=fileparts(qpdfiles{nn});

    pd_raw = load(qpdfiles{nn});     
    pd_data(nn).FullName = qpdfiles{nn};
    pd_data(nn).Params = params(nn);
    pd_data(nn).Data = pd_raw;
    disp(b);
end

%%
assignin('base','pd_data',pd_data);
assignin('base','pd_opts',pd_opts);
pd_analysis
end

