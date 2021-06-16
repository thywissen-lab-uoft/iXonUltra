function grabMagnetometer
%flowRatePlotter Summary of this function goes here
%   Detailed explanation goes here

% Question : How much to "GUIify" this interface?
% To do:
% File parser over many CSV files
% GUI / command line options for different variables
% GUI / command for a range of times
% Choosing y limits

% Root directory of logs
fldr='Y:\LabJack\Logging\Magnetometer';

% Start and End times to plot YYYY mm DD HH MM SS
t1=[2021 06 16 15 20 30];
t2=[2021 06 16 16 10 0];

% Resampling time in minutes (for average)
dt=0;

flim=[0 2];

% Load the logs
tic
rawTbl=loadLogs(t1,t2);
toc



% Resample the data using average to make plotting easier over long times
if dt
    tic
    pTable=retime(rawTbl,'regular','mean','timestep',minutes(dt));
    toc
else
    pTable=rawTbl;
end

hF=figure;
hF.Position=[50 100 1600 400];
set(hF,'color','w');
fnames=pTable.Properties.VariableNames;


C=500; % mG/V (1V/50uT according to spec sheet)

for kk=1:length(fnames)
     ps(kk)=plot(pTable.Time,C*pTable.(fnames{kk}),'linewidth',2);
     hold on
end

ylabel('field (mG)');
set(gca,'fontsize',12);

legend({'500 mG/V'});

% ylim([-85 -40]);
% ylim([-100 -40]);

xlim([datetime(t1) datetime(t2)]);
% 
function out=loadLogs(t1,t2)
    t1=datenum(t1);
    t2=datenum(t2);      
    out=readLog(makeFileName(t1));
    t1=t1+1;    
    while floor(t1)<=floor(t2)
%         disp(t1)
        str=makeFileName(t1);
        if exist(str,'file')      
            thisdata=readLog(str);
            if ~isempty(thisdata)
                out=[out;thisdata];
            end
        end
        t1=t1+1;          
    end 
end

function str=makeFileName(t)
   tV=datevec(t);       
   str=[fldr filesep num2str(tV(1)) filesep num2str(tV(1)) '.' num2str(tV(2),'%02.f') filesep num2str(tV(2),'%02.f') '_' datestr(t,'dd') '.csv'];
end

end

function out=readLog(fname)
% Could use readtable, but it is slightly slower than text scan.
% Over many csv files with will add up (could consider going to SQL to
% further reduce time?)
    
    fprintf('Reading file...')
    
    if exist(fname)    
        fid=fopen(fname);
        disp(fname)
        hdr=textscan(fgetl(fid),'%s','delimiter',',');
        hdr=hdr{1};nhdr=length(hdr);
        fmt=['%q',repmat('%f',1,nhdr-1)];
        data=textscan(fid,fmt,'delimiter',',');
        fclose(fid);  

        fprintf('Converting string to date...');
        data{:,1}=datetime(data{:,1},'InputFormat','MM/dd/yyyy, HH:mm:ss');    
        disp('done');


        fprintf('Making time table object ... ');
        out=timetable(data{:,1},data{:,2:end},'VariableNames',hdr(2:end));
        disp('done');    
    else
        disp('no file');
        out=[];
    end

end

