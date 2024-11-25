clear composite_data
composite_data = struct;
index=1;

%% Saving and Output

output_filename = '2024_11_shake_2.5Er_190G';
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%% Define the Runs

%vary field 50 ms mod ramp 11/11 evap to 65 mW
composite_data(index).Name = '11/11-12';
composite_data(index).Description = '2.5Er, 50 ms mod ramp, 0.4 V, 190 G, 65mW';
composite_data(index).Runs= [
        2024 11 11 10;
        2024 11 11 11;
        2024 11 11 12;
        2024 11 11 13;
        2024 11 11 14;
        2024 11 12 01;
        2024 11 12 02;
        2024 11 12 03;
        2024 11 12 04;
        2024 11 12 05;
        2024 11 12 06;
       
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