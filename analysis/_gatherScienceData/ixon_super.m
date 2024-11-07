function ixon_super(composite_data,composite_opts)

if nargin==1
    composite_opts=struct;
end

if ~isfield(composite_opts,'do_ixon_main')
    composite_opts.do_ixon_main = 1;
end

if ~isfield(composite_opts,'do_ixon_bin_analysis')
    composite_opts.do_ixon_bin_analysis = 1;
end

if ~isfield(composite_opts,'do_ixon_dig_analysis')
    composite_opts.do_ixon_dig_analysis = 1;
end
try
for super_composite_data_index=1:length(composite_data)
    dir_list = ixon_findRunDirectory(composite_data(super_composite_data_index).Runs);
    
    for super_directory_index=1:length(dir_list)
        ixon_auto_dir = 0;
        imgdir = dir_list{super_directory_index};
        if composite_opts.do_ixon_main
            ixon_main;
        end

        bin_auto_file = 0;
        imgdir = dir_list{super_directory_index};
        if composite_opts.do_ixon_bin_analysis
            filename = fullfile(dir_list{super_directory_index},'figures','bindata.mat');
            bin_imgdir = fullfile(dir_list{super_directory_index},'figures');
            ixon_bin_analysis;
        end

        dig_auto_file = 0;
        imgdir = dir_list{super_directory_index};
        if composite_opts.do_ixon_dig_analysis
            dig_imgdir = fullfile(dir_list{super_directory_index},'figures');
            filename = fullfile(dir_list{super_directory_index},'figures','digdata.mat');
            ixon_dig_analysis;
        end
    end   
end

catch ME

    keyboard

end



end

