   % y = 2024;

% YYYY=datestr(d,'YYYY');
% mm = datestr(d,'mm');
% dd = datestr(d,'dd');

    src='X:\LabJackLogs\ODTQPD';

YYYY = '2024';
mm = '11';

dd={'01','02','03','04','05','06','07','07','09','10','11','12'};

x=[];
y=[];
z=[];
t=datetime;
t(1)=[];

for kk=1:length(dd)    
    pd_dir = fullfile(src,YYYY,[YYYY '.' mm],[mm '.' dd{kk}]);
    flist=dir([pd_dir filesep 'ODTQPD_*.mat']);
    for nn=1:10:length(flist)
        disp(flist(nn).name);
        fullfilename = fullfile(pd_dir,flist(nn).name);
        data=load(fullfilename);
        x(end+1)=mean([data.data(end-100:end,7)]);
        y(end+1)=mean([data.data(end-100:end,8)]);
        z(end+1)=mean([data.data(end-100:end,9)]);
        t(end+1)=datetime(datenum(data.AcquisitionTime,'YYYY-mm-dd_HH-MM-SS'),'convertfrom','datenum');    end
end

% pd_dir = fullfile(src,YYYY,[YYYY '.' mm],[mm '.' dd]);
% flist=dir([pd_dir filesep 'ODTQPD_*.mat']);
% filename=flist(end).name;
% 


%%

t(1);

figure;
a1=subplot(311);
plot(t,x,'.');
ylim([-.5 6])
ylabel('x latt')

a2=subplot(312);
plot(t,y,'.');
ylim([-.5 9])
ylabel('y latt')

a3=subplot(313);
plot(t,z,'.');
ylim([-.5 6])
ylabel('z latt')
linkaxes([a1 a2 a3],'x')
