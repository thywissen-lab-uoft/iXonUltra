function [ixondata,pd_summary] = AnalyzeIxonQPD(ixondata,saveDir,opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matches QPD files to a list of iXon data of the form ixondata.digdata
% Analyzes each run to get average powers, fit parameters for data list

% opts fields are:
% doPlot (plots QPD analyis)
% doSave (saves plots to iXon folder)
% isMac (sets source dir.)
% isRemote 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Display

disp('Analyzing QPD traces...')

%% Initialize options

if nargin~=3
    opts.doPlot = 1;
    opts.isMac = 0;
    opts.doSave = 1;
    opts.isRemote = 0;
end

%% Data Root
% Find the data source file

if opts.isMac       
    pdsrc =  '/Volumes/main/LabJackLogs/ODTQPD'; % Mac
else
    pdsrc = 'X:\LabJackLogs\ODTQPD'; %Windows
end


%% Load the QPD files
qpdfiles=strings(0);

for jj=1:length([ixondata.Params])

    % Date string for run
    fdate = datestr(ixondata(jj).Params.ExecutionDateStr,'YYYY-mm-dd_HH-MM-SS');
    
    % Directory of QPD file
    pddir = fullfile(pdsrc,fdate(1:4),[fdate(1:4) '.' fdate(6:7)],[fdate(6:7) '.' fdate(9:10)]);
    
    % QPD filename
    fname = ['ODTQPD_' fdate '.mat'];
    qpdfiles(end+1) = fullfile(pddir,fname);

end

%% Analyze the QPD traces

clear qpd_data
qpd_data = photodiode_analyze(qpdfiles);    % Analyze a single trace

%% Assign single QPD analysis trace to ixondata
for nn=1:length(qpdfiles)
    ixondata(nn).qpd_data = qpd_data(nn);           % Assign output of single trace to ixondata
end

%% Perform the analysis on the total qpd data

[pd_summary,f_trace,f_mod,f_sum] = photodiode_collect_analysis(qpd_data,opts);


if opts.doSave
    saveas(f_trace,[saveDir filesep 'qpd_traces.png']);
    saveas(f_mod,[saveDir filesep 'qpd_modulation.png']);
    saveas(f_sum,[saveDir filesep 'qpd_power.png']);
    filename=fullfile(saveDir,'pd_summary.mat');
    save(filename,'-struct','pd_summary');    
end



end