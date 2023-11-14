function ixon_ps_qp_compare(ixondata)

src='X:\LabJackLogs\CATS';
Iall = zeros(length(ixondata),2);
for kk=1:length(ixondata)
    t = ixondata(kk).Params.ExecutionDate;
   
    yyyy= datestr(t,'yyyy');
    mm= datestr(t,'mm');
    dd= datestr(t,'dd');
    HH= datestr(t,'HH');
    MM= datestr(t,'MM');
    SS= datestr(t,'SS');

%     yyyy = f(1:4);
%     mm=f(6:7);
%     dd=f(9:10);
%     HH=f(12:13);
%     MM = f(15:16);
%     SS = f(18:19);
%     fdata=f;
    tdata=[str2num(yyyy) str2num(mm) str2num(dd) str2num(HH) str2num(MM) str2num(SS)];
    
    mydir = fullfile(src,yyyy,[yyyy '.' mm],[mm '.' dd]);    
    names = dir(fullfile(mydir,'*.mat'));
    names={names.name};
    names=sort(names);
    

    
    for jj = 1:length(names)
        if contains(names{jj},'extra')
           continue 
        end
        f=names{jj};
        f=f(6:end);        
        yyyy = f(1:4);
        mm=f(6:7);
        dd=f(9:10);
        HH=f(12:13);
        MM = f(15:16);
        SS = f(18:19);
        
        tcats=[str2num(yyyy) str2num(mm) str2num(dd) str2num(HH) str2num(MM) str2num(SS)];
        dT=(datenum(tcats)-datenum(t))*24*60;
        if abs(dT)<3            
            str = fullfile(src,yyyy,[yyyy '.' mm],[mm '.' dd],names{jj});    

            
            cats = load(str);
            
            t1 = cats.SequencerTime;
            
            
            t2 = datestr(t,'yyyy-mm-dd_HH-MM-SS');
            
            if isequal(t1,t2)
                break
            end
            
        end      

    end
    
    x = cats.t;
    y = cats.data(:,end-1);
    y0 = mean(y(end-500:end));
    y= y-y0;
    y = y*1000/50;
    
    ta = 59.18;
    tb = 59.8;
    
    ii = logical([x>ta].*[x<tb]);
    
    Ithis = mean(y(ii));
    dIthis = std(y(ii));
%     keyboard
    Iall(kk,1)=Ithis;
        Iall(kk,2)=dIthis;

end
P=[ixondata.Params];
itot = [P.qgm_plane_IQP];
keyboard

end

