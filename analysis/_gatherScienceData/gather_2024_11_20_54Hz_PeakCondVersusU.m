clear composite_data
composite_data = struct;
index=1;

%% 2024/11/20-2024/11/21
% 50 ms mod ramp 2024/11/20-2024/11/21 vary field
composite_data(index).Name = '11/20 0.4V 54 Hz';
composite_data(index).Description = '2.5Er, 54 Hz, 0.4V, 65 mW, vary field';
composite_data(index).Type = 'peak';
composite_data(index).Runs =[     
    2024 11 20 05;
    2024 11 20 06;
    2024 11 21 01;
    2024 11 21 02;
    2024 11 21 03;
    2024 11 21 04;
    2024 11 21 05;
    2024 11 21 06;
    2024 11 21 07;
    2024 11 21 08;
    2024 11 21 09;
    ];
index=index+1;

%% Redo Analysis
do_redo_analysis = 1;    % Do you want to run analysis on it?

if do_redo_analysis
    opts=struct;
    opts.do_ixon_main           = 0;   % ixon_main
    opts.do_ixon_bin_analysis   = 0;   % ixon_bing
    opts.do_ixon_dig_analysis   = 1;   % ixon_dig
    ixon_super(composite_data,opts)
end

%% Gather Data
composite_data = gatherCompositeData(composite_data);
%% Upload

doUpload = true;

GDrive_root =['G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name = '2024_11_20 Peak Cond 54 Hz Versus U';
saveDir = fullfile(GDrive_root,output_folder_name);

if doUpload
    try
        if ~exist(GDrive_root,'dir');mkdir(GDrive_root);end
        if ~exist(saveDir,'dir');mkdir(saveDir);end
         gFile = fullfile(saveDir,'composite_data.mat');

        disp(gFile);
        fprintf('upload to google drive ...');        
        save(gFile,'composite_data'); 
        disp('done!')   
    catch ME
        error('OH NO I COUD NOT UPL<OD');
    end
end 

