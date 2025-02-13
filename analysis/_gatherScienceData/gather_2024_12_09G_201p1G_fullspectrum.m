clear composite_data
composite_data = struct;
index=1;
%% Introduction
%% 12/09-12/10

composite_data(index).Name = '2024_12_09 201.1 G 66.5 mW 5 Er pulse spectrum';
composite_data(index).Description = '2.5 Er, 50 ms mod ramp, 66.5 mW 5 Er pulse spectrum';
% composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [    
    2024 12 09 02;
    2024 12 09 03;
    2024 12 09 04;
    2024 12 09 05;
    2024 12 09 06;
    2024 12 10 01;
    2024 12 10 02;
    2024 12 10 03;
    2024 12 10 04;
    2024 12 10 05;
    2024 12 10 06;
    2024 12 10 07;
    2024 12 10 08;
    2024 12 10 09;
    2024 12 10 10;
    2024 12 10 11;
    2024 12 10 12;
    2024 12 10 13;
    2024 12 10 14;
    ];

index=index+1;
%% Redo Analysis
do_redo_analysis = 0;    % Do you want to run analysis on it?

if do_redo_analysis
    opts=struct;
    opts.do_ixon_main           = 1;   % ixon_main
    opts.do_ixon_bin_analysis   = 1;   % ixon_bing
    opts.do_ixon_dig_analysis   = 1;   % ixon_dig
    ixon_super(composite_data,opts)
end

%% Gather Data
composite_data = gatherCompositeData(composite_data);

%% Upload

doUpload = true;

GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name = '2024_12_09 201.1 G 66.5 mW 5 Er pulse spectrum';
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