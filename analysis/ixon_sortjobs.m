mydir = 'X:\Data\2024\2024.03\03.04\03 overnight, 199.4 G, 57mW AC shake';
names = dir([mydir filesep '*.mat']);
names = {names.name};

daydir = 'X:\Data\2024\2024.03\03.04';
% dir(daydir


dir1 =  'X:\Data\2024\2024.03\03.04\03a tiltcheck';
dir2 = mydir;
if ~exist(dir1,'dir')
    mkdir(dir1);
end
if ~exist(dir2,'dir')
    mkdir(dir2);
end
for kk=1:length(names)
   filename = fullfile(mydir,names{kk});
   d = load(filename);
    index = d.data.Flags.plane_selection.dotilt;
    
    
    index = [d.data.Params.xdt_evap1_power > .07];
    if index
        movefile(filename,dir1);
    else        
        amp = d.data.Params.conductivity_ODT2_mod_amp;
        freq = d.data.Params.conductivity_mod_freq;  
        pow = d.data.Params.xdt_evap1_power*1000;
        field = d.data.Params.conductivity_FB_field;
        u=d.data.Params.latt_depth_load(1);
        rt = d.data.Params.conductivity_mod_ramp_time;
        
        dirname = ['shake plane ' num2str(freq) ' Hz, ' ...
            num2str(field) ' G, ' ...
            num2str(u) ' Er, ' ...
            num2str(pow) ' mW, ' ...
            num2str(amp,'%.1f') ' V ODT2,' ...
            num2str(rt) ' ms ramp time'];        
        dirnametot = [daydir filesep dirname];
        if ~exist(dirnametot,'dir')
           mkdir(dirnametot); 
        end
        movefile(filename,dirnametot);
    end
end