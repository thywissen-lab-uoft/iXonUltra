clear composite_data
composite_data = struct;
index=1;

%% Saving and Output

output_filename = '2024_11_shake_2.5Er_201.1G';
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%% Define the Runs


% % 201.1 G (missing QPD analysis)
% composite_data(index).Name = '201.1 G 03/14';
% composite_data(index).Description = 'uh';
% composite_data(index).Runs = [
%     2024 03 14 06;
%     2024 03 14 07;
%     2024 03 14 08;
%     2024 03 14 09;
%     2024 03 14 10;
%     2024 03 14 11;
%     2024 03 14 12;
%     2024 03 14 14;
%     2024 03 14 15;
%     2024 03 14 16;
%     2024 03 14 17;
%     2024 03 14 18;
%     2024 03 14 19;
%     2024 03 14 20;
%     2024 03 14 21;
%     2024 03 14 22;
%     ];
% index=index+1;

% 201.1 G high field evap 50 ms mod ramp 10/24 lower drive evap to 65 mW
composite_data(index).Name = '201.1 G 10/24';
composite_data(index).Description = '201.1 G high field evap 50 ms mod ramp 10/24 lower drive evap to 65 mW';
composite_data(index).Runs=  [ 
    2024 10 24 07;
    2024 10 24 08;
    2024 10 24 09;
    2024 10 24 10;
    2024 10 24 11;
    2024 10 24 12;
    2024 10 24 13;
    2024 10 24 14;
    2024 10 24 15;
    2024 10 24 16;
    2024 10 24 17;
    2024 10 24 18;
    2024 10 24 19;
    ];
index=index+1;

%201.1 G high field evap 50 ms mod ramp 10/31 evap to 65.5 mW, hold at
%201.1G for 200 ms before shaking
composite_data(index).Name = '201.1 G 10/31';
composite_data(index).Description = ['201.1 G high field evap 50 ms mod ramp 10/31 evap to 65.5 mW ' ...
    ', hold at 201.1G for 200 ms before shaking'];
composite_data(index).Runs= [ 
    2024 10 31 07;
    2024 10 31 09;
    2024 10 31 11;
    2024 10 31 13;
    2024 10 31 15;
    2024 10 31 17;
    2024 10 31 19;
    2024 10 31 21;
    2024 10 31 23;
    2024 10 31 25;
    ];
index=index+1;

%201.1 G high field evap 50 ms mod ramp 11/01 evap to 64 mW
composite_data(index).Name = '201.1 G 11/01';
composite_data(index).Description = '201.1 G high field evap 50 ms mod ramp 11/01 evap to 64 mW';
composite_data(index).Runs= [ 
    2024 11 01 07;
    2024 11 01 08;
    2024 11 01 09;
    2024 11 01 10;
    2024 11 01 11;
    2024 11 01 12;
    2024 11 01 13;
    2024 11 01 14;
    2024 11 01 15;
    2024 11 01 16;
    2024 11 01 17;
    2024 11 01 18;
    ];
index=index+1;

%201.1 G high field 50 ms mod ramp 11/05 evap to 65.5 mW
composite_data(index).Name = '201.1 G 11/05';
composite_data(index).Description = '201.1 G high field 50 ms mod ramp 11/05 evap to 65.5 mW';
composite_data(index).Runs= [ 
    2024 11 05 09;
    2024 11 05 10;
    2024 11 05 11;
    2024 11 05 12;
    2024 11 05 13;
    2024 11 05 14;
    2024 11 05 15;
    2024 11 05 16;
    2024 11 05 17;
    2024 11 06 01;
    2024 11 06 02;
    2024 11 06 03;
    2024 11 06 04;
    2024 11 06 05;
    2024 11 06 06;
    2024 11 06 07;
    2024 11 06 08;
    2024 11 06 09;
    2024 11 06 10;
    2024 11 06 11;
    2024 11 06 12;
    2024 11 06 13;
    2024 11 06 14;
    2024 11 06 15;
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