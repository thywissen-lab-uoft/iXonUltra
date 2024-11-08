clear composite_data
composite_data = struct;
index=1;

%% Saving and Output

output_filename = '2024_11_peakCondVTemp_2.5Er_201.1G';
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%% Define the Runs

%201.1 G high field 50 ms mod ramp 11/05 evap to 65.5 mW
composite_data(index).Name = '201.1 G 11/07 11/08';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 51 Hz 0.9V, vary lattice pulse depth ';
composite_data(index).Runs= [ 
    2024 11 07 18;
    2024 11 07 19;
    2024 11 07 20;
    2024 11 07 21;
    2024 11 07 22;
    2024 11 07 23;
    2024 11 07 24;
    2024 11 07 25;
    2024 11 07 26;
    2024 11 07 27;
    2024 11 07 28;
    2024 11 07 29;
    2024 11 08 01;
    2024 11 08 02;
    2024 11 08 03;
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