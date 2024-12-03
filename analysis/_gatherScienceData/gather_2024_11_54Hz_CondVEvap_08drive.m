clear composite_data
composite_data = struct;
index=1;

%% Saving and Output

output_filename = '2024_11_peakCondVEvap_2.5Er_201.1G_0.8Vdrive';
doUpload = 1;
GDrive_root =['G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%% Define the Runs

%201.1 G high field 50 ms mod ramp 11/25-11/26 vary evap depth
composite_data(index).Name = '201.1 G 11/25 26';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 54 Hz 0.8V, vary evap depth ';
composite_data(index).Type = 'peak';
composite_data(index).Runs= [ 
        2024 11 25 22; 
        2024 11 25 23; % questionable
        2024 11 25 24; % bad SNR
        2024 11 26 01; 
        2024 11 26 02; 
        2024 11 26 03; % good 68.5 mW
        2024 11 26 04; % atom number varies
        2024 11 26 05;
        2024 11 26 06; % questionable fit - below 0.5um, 69.5 mW
        2024 11 26 07; % questionable fit
        2024 11 26 08; % questionable fit 
        2024 11 26 09;
        2024 11 26 10; % questionable fit 70 mW

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
