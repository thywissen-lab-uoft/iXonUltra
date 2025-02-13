function qpd_super(composite_data,composite_opts)


for super_composite_data_index=1:length(composite_data)
    dir_list = ixon_findRunDirectory(composite_data(super_composite_data_index).Runs);
    
    for super_directory_index=1:length(dir_list)
        
        try 
            filename = fullfile(dir_list{super_directory_index},'figures','digdata.mat');
            bd = load(filename);    
        catch
            filename = fullfile(dir_list{super_directory_index},'figures','ixon_boxdata.mat');
            bd = load(filename);
            bd = bd.ixon_boxdata;
        end
        
        D = [bd.Params.ExecutionDate];
        composite_opts.X = [bd.X];
        composite_opts.xVar = bd.xVar;
        composite_opts.Frequency = unique([bd.Params.conductivity_mod_freq]);
        composite_opts.RampTime  = unique([bd.Params.conductivity_mod_ramp_time]);
        
        composite_opts.saveDir =  fullfile(dir_list{super_directory_index},'figures');

        try
            close all;
            [figs,output] = qpd_main(D,composite_opts) ;
            
            for kk=1:length(figs)
                  if ~isempty(figs{kk})
                    ixon_saveFigure2(figs{kk},...
                     figs{kk}.Name,composite_opts);  
                  end
            end
          
            try if ~exist(composite_opts.saveDir,'dir');mkdir(composite_opts.saveDir);end;end
            
            filename = fullfile(composite_opts.saveDir,'qpd.mat');
            disp(['Saving ' filename ' ...']);
            save(filename, '-struct','output');      

        catch ME
           warning(getReport(ME),'extended','hyperlinks','on');
           disp('DID YOU MAKE SURE TO ADD THE AUXLIARY_ANALYSIS GIT REPOSITORY IS LOADED?');
        end
        
    end
    
end

end

