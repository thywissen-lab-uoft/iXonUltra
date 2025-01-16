clear composite_data
composite_data = struct;
index=1;
%% Introduction
%% 2024/12/05

 % 2.5 Er, 50 ms mod ramp, 200G, 64-70 mW, 4 Er pulse
composite_data(index).Name = '2024_12_05 201.1 G 64.5 mW spectrum';
composite_data(index).Description = '2.5 Er, 50 ms mod ramp, 201.1 G, 64.5 mW';
% composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [    
    2024 12 05 01;
    2024 12 05 02;
    2024 12 05 03;
    2024 12 05 04;
    2024 12 05 05;
    2024 12 06 01;
    2024 12 06 02;
    2024 12 06 03;
    2024 12 06 04;
    2024 12 06 05;
    2024 12 06 06;
    2024 12 06 07;
    2024 12 06 08;
    2024 12 06 09;
    2024 12 06 10;
    2024 12 06 11;
    2024 12 06 12;
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


saveDir =['G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_name = '2024_12_05_full_spectrum_190G_T1.mat';


if doUpload
    try
        if ~exist(saveDir,'dir');mkdir(saveDir);end
         gFile = fullfile(saveDir,output_name);

        disp(gFile);
        fprintf('upload to google drive ...');        
        save(gFile,'composite_data'); 
        disp('done!')   
    catch ME
        error('OH NO I COUD NOT UPL<OD');
    end
end 