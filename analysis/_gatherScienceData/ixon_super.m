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

for kk=1:length(composite_data)
    dir_list = ixon_findRunDirectory(composite_data(kk).Runs);
    
    for nn=1:length(dir_list)
        ixon_auto_dir = 0;
        bin_auto_file = 0;
        dig_auto_file = 0;
        imgdir = dir_list{nn};

        if composite_opts.do_ixon_main
            imgdir = dir_list{nn};
            ixon_main;
        end

        if composite_opts.do_ixon_bin_analysis
            filename = fullfile(dir_list{nn},'figures','bindata.mat');
            bin_imgdir = fullfile(dir_list{nn},'figures');
            ixon_bin_analysis;
        end

        if composite_opts.do_ixon_dig_analysis
            dig_imgdir = fullfile(dir_list{nn},'figures');
            filename = fullfile(dir_list{nn},'figures','digdata.mat');
            ixon_dig_analysis;
        end
    end   
end


end

