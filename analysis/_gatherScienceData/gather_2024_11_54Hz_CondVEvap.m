clear composite_data
composite_data = struct;
index=1;

%% Saving and Output

output_filename = '2024_11_peakCondVEvap_2.5Er_201.1G';
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%% Define the Runs

%201.1 G high field 50 ms mod ramp 11/22-11/23 vary evap depth
composite_data(index).Name = '201.1 G 11/22 23';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 54 Hz 0.4V, vary evap depth ';
composite_data(index).Runs= [ 
%         2024 11 22 20; % bad SNR
%         2024 11 22 21; % bad SNR, number fluctuating
%         2024 11 22 22; % bad SNR
        2024 11 22 23;
        2024 11 22 24;
        2024 11 22 25;
        2024 11 22 26;
        2024 11 22 27; % questionable
%         2024 11 22 28; % bad SNR
        2024 11 23 01;
        2024 11 23 02; % questionable
        2024 11 23 03;
        2024 11 23 04;
        
%         2024 11 23 05; % bad SNR
%         2024 11 23 06; % bad SNR
%         2024 11 23 07; % bad SNR
        2024 11 23 08; % questionable
        2024 11 23 09;
        2024 11 23 10;
        2024 11 23 11;
        2024 11 23 12;
        2024 11 23 13; % questionable, drifting atom number
        2024 11 23 14;
%         2024 11 23 15; % bad snr
        2024 11 23 16;
        2024 11 23 17;
%         2024 11 23 18; % bad SNR
        2024 11 23 19; % questionable
    ];
index=index+1;

%% Define the Runs

%201.1 G high field 50 ms mod ramp 11/22-11/23 vary evap depth
composite_data(index).Name = '201.1 G 11/22 23';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 54 Hz 0.4V, vary evap depth ';
composite_data(index).Runs= [ 
%         2024 11 22 20; % bad SNR
%         2024 11 22 21; % bad SNR, number fluctuating
%         2024 11 22 22; % bad SNR
        2024 11 22 23;
        2024 11 22 24;
        2024 11 22 25;
        2024 11 22 26;
        2024 11 22 27; % questionable
%         2024 11 22 28; % bad SNR
        2024 11 23 01;
        2024 11 23 02; % questionable
        2024 11 23 03;
        2024 11 23 04;
        
%         2024 11 23 05; % bad SNR
%         2024 11 23 06; % bad SNR
%         2024 11 23 07; % bad SNR
        2024 11 23 08; % questionable
        2024 11 23 09;
        2024 11 23 10;
        2024 11 23 11;
        2024 11 23 12;
        2024 11 23 13; % questionable, drifting atom number
        2024 11 23 14;
%         2024 11 23 15; % bad snr
        2024 11 23 16;
        2024 11 23 17;
%         2024 11 23 18; % bad SNR
        2024 11 23 19; % questionable
    ];
index=index+1;
%%
runs = [     
    2024 11 25 22;
    2024 11 25 23;
    2024 11 25 24;
    2024 11 26 01;
    2024 11 26 02;
    2024 11 26 03;
    2024 11 26 04;
    2024 11 26 05;
    2024 11 26 06;
    2024 11 26 07;
    2024 11 26 08;
    2024 11 26 09;
    2024 11 26 10;
    ];


composite_data(index).Name = '54 HZ 201.1 G 0.8V vary evap depth';
composite_data(index).Runs = runs;
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 54 Hz 0.8V, vary evap depth ';

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