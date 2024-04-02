
% Source directory of all data
srcdir = 'X:\Data\2024\2024.03\03.28\04 overnight 190 G and 195 G';


% Parent directory of data to sort
daydir = fileparts(srcdir); 


%% Find next available directory index
makeNewTiltDir = 1;
d = dir(daydir);
dfolders = d([d(:).isdir]);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));

ind = 0;
for kk = 1:length(dfolders)
    try
       nn = dfolders(kk).name ;
       bb = str2double(nn(1:2));
       if ~isnan(bb) && ~isinf(bb)
          ind = max([bb ind]); 
       end
       
       if isequal(strrep(nn(3:end),' ',''),'tiltcheck')
           makeNewTiltDir = 0;
           tiltdir = fullfile(daydir,nn);
       end
    end
end

%% Make tilt check folder

if makeNewTiltDir 
    dirtilt = fullfile(daydir,[num2str(ind,'%02.0f') ' ' 'tiltcheck']);
    if ~exist(dirtilt)
       mkdir(dirtilt); 
    end
end

%% Sort all files

% Get all file names
names = dir([srcdir filesep '*.mat']);
names = {names.name};

for kk=1:length(names)
    fprintf([num2str(kk) '/' num2str(length(names)) ' ']);
    fprintf(names{kk})
    fprintf('--> ');
    filename = fullfile(srcdir,names{kk});
    d = load(filename);
    data= d.data;
    
    P = data.Params;
    F = data.Flags;
    
    dotilt = F.plane_selection.dotilt;
    pow = 1e3*P.xdt_evap1_power;
    
    % ODT mode (0: off, 1: shake);
    m1 = F.conductivity_ODT1_mode;
    m2 = F.conductivity_ODT2_mode;    
    
    if dotilt
        [~,dirname,~]=fileparts(dirtilt);
        disp(dirname);
        movefile(filename,dirtilt);
    else
        if m1
            A2  = P.conductivity_ODT2_mod_amp;
            f   = P.conductivity_mod_freq;  
            B   = P.conductivity_FB_field;
            U   = P.latt_depth_load(1);
            tau = P.conductivity_mod_ramp_time; 
            dirname = ['shake plane ' num2str(f) ' Hz, ' ...
                num2str(round(B,2)) ' G, ' ...
                num2str(round(U,2)) ' Er, ' ...
                num2str(round(pow,2)) ' mW, ' ...
                num2str(A2,'%.1f') ' V ODT2, ' ...
                num2str(tau) ' ms ramp time'];   
        else 
            B   = P.conductivity_FB_field;
            U   = P.latt_depth_load(1);
            dirname = ['single plane ' num2str(round(B,2)) ...
               ' G,' num2str(round(U,2)) ' Er'];
        end
    
        fulldir = fullfile(daydir,[num2str(ind,'%02.0f') ' ' dirname]);
        if ~exist(fulldir)
            ind = ind+1;
            fulldir = fullfile(daydir,[num2str(ind,'%02.0f') ' ' dirname]);
            mkdir(fulldir);
        end        
        
        disp([num2str(ind,'%02.0f') ' ' dirname])
         movefile(filename,fulldir);

    end
    
    
    
%     disp(fulldir);
    

    
    % mod parameters
%     amp = data
    
%     fprintf(num2str(dotilt));
    
%     P = [d.

%     
%     dotilt = [d.data.Params.xdt_evap1_power > .07];
%     if dotilt
%         movefile(filename,dir1);
%     else        
%         amp = d.data.Params.conductivity_ODT2_mod_amp;
%         freq = d.data.Params.conductivity_mod_freq;  
%         pow = d.data.Params.xdt_evap1_power*1000;
%         field = d.data.Params.conductivity_FB_field;
%         u=d.data.Params.latt_depth_load(1);
%         rt = d.data.Params.conductivity_mod_ramp_time;
%         
%         dirname = ['shake plane ' num2str(freq) ' Hz, ' ...
%             num2str(field) ' G, ' ...
%             num2str(u) ' Er, ' ...
%             num2str(pow) ' mW, ' ...
%             num2str(amp,'%.1f') ' V ODT2,' ...
%             num2str(rt) ' ms ramp time'];        
%         dirnametot = [daydir filesep dirname];
%         if ~exist(dirnametot,'dir')
%            mkdir(dirnametot); 
%         end
%         movefile(filename,dirnametot);
%     end
%     disp(' ');
end