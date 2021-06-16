function [T,B]=grabMagnetometer(t1,t2)
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
% t1=[2021 06 16 15 20 30];
% t2=[2021 06 16 16 10 0];

C=500; % mG/V (1V/50uT according to spec sheet)

% Load the logs

tic
rawTbl=loadLogs(t1(1:6),t2(1:6));
toc

T=rawTbl.Time;
fname=rawTbl.Properties.VariableNames{1};
B=C*rawTbl.(fname);

i1=T>datetime(t1);
T=T(i1);
B=B(i1);

i2=T<datetime(t2);
T=T(i2);
B=B(i2);

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

