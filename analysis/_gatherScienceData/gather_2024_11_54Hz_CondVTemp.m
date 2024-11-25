clear composite_data
composite_data = struct;
index=1;

%% Saving and Output

output_filename = '2024_11_54Hz_CondVTemp_2.5Er_201.1G';
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%% Define the Runs

%201.1 G high field 50 ms mod ramp 11/05 evap to 65.5 mW
composite_data(index).Name = '201.1 G 11/21-22';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 54 Hz 0.4V, vary lattice pulse depth ';
composite_data(index).Runs= [ 
        2024 11 21 10;
%         2024 11 21 11; % bad cond SNR
        2024 11 21 12;
        2024 11 21 13;
        2024 11 21 14;
        2024 11 21 15;
        2024 11 21 16;
%         2024 11 21 17; % bad cond SNR
        2024 11 21 18;
%         2024 11 22 01;
        2024 11 22 02;
        2024 11 22 03;
        2024 11 22 04;
%         2024 11 22 05;
        2024 11 22 06;
        2024 11 22 07;
        2024 11 22 08;
        2024 11 22 09;
        2024 11 22 10;
%         2024 11 22 11;
        2024 11 22 12;
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