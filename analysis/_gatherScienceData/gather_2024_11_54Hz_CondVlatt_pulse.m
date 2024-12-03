clear composite_data
composite_data = struct;
index=1;

%% Saving and Output

output_filename = '2024_11_peakCondVEvap_2.5Er_201.1G';
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%% Define the Runs

%201.1 G high field 50 ms mod ramp 11/22-11/23 vary evap depth
composite_data(index).Name = '201.1 G 11/23 24';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 54 Hz 0.4V, vary evap depth, pulse 4 ER lattice ';
composite_data(index).Runs= [ 
        2024 11 23 20; % atom number fluct affecting fit, questionalbe SNR
        2024 11 23 21; % atom number drifting 
        2024 11 23 22; % bad SNR
        2024 11 23 23; % bad SNR
        2024 11 23 24; % bad SNR
        2024 11 23 25; % bad SNR
        2024 11 23 26; % bad SNR
        2024 11 23 27;
        2024 11 23 28;
        2024 11 24 01;
        2024 11 24 02; % questionable fit
        2024 11 24 03; % bad SNR
        2024 11 24 04;
   
        
        2024 11 24 05;
        2024 11 24 06;
        2024 11 24 07;
        2024 11 24 08;
        2024 11 24 09;
        2024 11 24 10;
        2024 11 24 11;
        2024 11 24 12;
        2024 11 24 13;
        2024 11 24 14;
        2024 11 24 15;
        2024 11 24 16;
        2024 11 24 17;
    ];
index=index+1;

%% Redo Analysis
do_redo_analysis = 0;    % Do you want to run analysis on it?

if do_redo_analysis
    opts=struct;
    opts.do_ixon_main           = 0;   % ixon_main
    opts.do_ixon_bin_analysis   = 0;   % ixon_bing
    opts.do_ixon_dig_analysis   = 1;   % ixon_dig
    ixon_super(composite_data,opts)
end

%% Gather Data
composite_data = gatherCompositeData(composite_data);

%% Upload the data

 try
     if ~exist(GDrive_root,'dir')
     mkdir(GDrive_root);
     end
 end
 
 if  doUpload && exist(GDrive_root,'dir')   
     fprintf('upload to google drive ...');
    gFile = fullfile(GDrive_root,[output_filename '.mat']);
    save(gFile,'composite_data'); 
    disp('done!')
 end