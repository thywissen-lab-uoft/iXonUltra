

%% Saving and Output
clear composite_data
composite_data = struct;
index=1;
output_filename = '2024_11_shake_2.5Er_201.1G_single_tone';
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%201.1 G high field 50 ms mod ramp 11/05 evap to 65.5 mW
composite_data(index).Name = '201.1 G 11/06';
composite_data(index).Runs= [ 
       2024 11 06 18;
       2024 11 06 19;
       2024 11 06 20;
       2024 11 07 01;
    2024 11 07 02;


    ];
index=index+1;
%% Saving and Output
clear composite_data
composite_data = struct;
index=1;
output_filename = '2024_11_shake_2.5Er_201.1G_single_tone_2';
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%201.1 G high field 50 ms mod ramp 11/05 evap to 65.5 mW
composite_data(index).Name = '201.1 G 11/06';
composite_data(index).Runs= [ 
    2024 11 07 03;
    2024 11 07 04;
    2024 11 07 05;
    2024 11 07 06;
    2024 11 07 07;
    2024 11 07 08;
    2024 11 07 09;
    2024 11 07 10;
    2024 11 07 11;

    ];
index=index+1;

%% Redo Analysis
do_redo_analysis = 1;    % Do you want to run analysis on it?

if do_redo_analysis
    opts=struct;
    opts.do_ixon_main           = 1;   % ixon_main
    opts.do_ixon_bin_analysis   = 1;   % ixon_bing
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