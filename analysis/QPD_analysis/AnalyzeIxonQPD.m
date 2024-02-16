function [ixondata,qpd_out] = AnalyzeIxonQPD(ixondata,saveDir,opts)

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

for nn=1:length(qpdfiles)
    ixondata(nn).qpd_data = qpd_data(nn);           % Assign output of single trace to ixondata
end

%% Perform the analysis on the qpd data

[hF,qpd_summary_data] = photodiode_collect_analysis(qpd_data,opts);


% if opts.doSave
%     saveas(fmsum,[saveDir filesep 'qpd_modulation_summary.png']);
%     saveas(fmsum_um,[saveDir filesep 'qpd_modulation_summary_um.png']);
%     saveas(fpsum,[saveDir filesep 'qpd_power_summary.png']);
% 
%     filename=fullfile(saveDir,'qpddata.mat');
%     save(filename,'qpd_out');    
% end



end