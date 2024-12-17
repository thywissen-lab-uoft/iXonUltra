clear composite_data
composite_data = struct;
index=1;

%% Saving and Output
output_filename = '2024_12_2.5ER_201.1G_full_spectrum_plane_shift';
doUpload = 1;
GDrive_root =['G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];

%% Define the Runs

%201.1 G high field 50 ms mod ramp 12/09-12/10 central plane full spectrum

composite_data(index).Name = '2024_12_09 201.1 G 66.5 mW 5 Er full spectrum';
composite_data(index).Description = '2024_12_09 201.1 G 66.5 mW 5 Er full spectrum, central plane N=-3';
composite_data(index).Type = 'spectrum';
composite_data(index).Runs= [ 
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
%%
%201.1 G high field 50 ms mod ramp 12/1-12/11 5 shift from central plane full spectrum

composite_data(index).Name = '2024_12_10 201.1 G 66.5 mW 5 Er full spectrum plane shift 5';
composite_data(index).Description = '2024_12_1 201.1 G 66.5 mW 5 Er full spectrum, 5 from central plane N=-5';
composite_data(index).Type = 'spectrum';
composite_data(index).Runs= [ 
        2024 12 10 15;
        2024 12 10 16;
        2024 12 10 17;
        % 2024 12 10 18; % between planes?
        2024 12 10 19;
        2024 12 10 20; 
        % 2024 12 10 21;% between planes
        % 2024 12 10 22; % between planes
        2024 12 10 23;
        2024 12 10 24;
    ];
index=index+1;

%%
%201.1 G high field 50 ms mod ramp 12/10-12/11 central plane full spectrum

composite_data(index).Name = '2024_12_1 201.1 G 66.5 mW 5 Er full spectrum';
composite_data(index).Description = '2024_12_1 201.1 G 66.5 mW 5 Er full spectrum, central plane N=-5';
composite_data(index).Type = 'spectrum';
composite_data(index).Runs= [ 
        2024 12 10 25;
        2024 12 10 26;
        2024 12 11 01;
        2024 12 11 02;
        2024 12 11 03;
        2024 12 11 04;
        2024 12 11 05;
        2024 12 11 06;
        2024 12 11 07;
        2024 12 11 08;
    ];
index=index+1;

%%

%%%%% ALLL BAD
% %201.1 G high field 50 ms mod ramp 12/10-12/11 5 shift from central plane full spectrum
% 
% composite_data(index).Name = '2024_12_1 201.1 G 66.5 mW 5 Er full spectrum';
% composite_data(index).Description = '2024_12_1 201.1 G 66.5 mW 5 Er full spectrum, central plane N=-7';
% composite_data(index).Type = 'spectrum';
% composite_data(index).Runs= [ 
%         2024 12 11 09; % out of focus
%         2024 12 11 10; % selecting between planes??
%         2024 12 11 11;
%         2024 12 11 12;
%         2024 12 11 13;
%         2024 12 11 14;
%         2024 12 11 15;
%         2024 12 11 16;
%         2024 12 11 17;
%         2024 12 11 18;
% 
% 
%     ];
% index=index+1;
% 

%% 12/11-12/12

% 2.5 Er, 50 ms mod ramp, 200G, 64-70 mW, 4 Er pulse
composite_data(index).Name = '2024_12_11 201.1 G 66.5 mW 5 Er pulse spectrum, -12 plane shift';
composite_data(index).Description = '2.5 Er, 50 ms mod ramp, 66.5 mW 5 Er pulse spectrum, -12 plane shift';
% composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [    
    2024 12 11 19;
    2024 12 11 20;
    2024 12 11 21;
    2024 12 11 22;
    2024 12 11 23;
    2024 12 11 24;
    2024 12 11 25;
    2024 12 12 01;
    2024 12 12 02;
    2024 12 12 03;
    ];

index=index+1;

% 2.5 Er, 50 ms mod ramp, 201.1 G, 66.5 mW, 5 Er pulse
composite_data(index).Name = '2024_12_11 201.1 G 66.5 mW 5 Er pulse spectrum, -8 plane shift';
composite_data(index).Description = '2.5 Er, 50 ms mod ramp, 66.5 mW 5 Er pulse spectrum, -8 plane shift';
% composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [    
    2024 12 12 04;
    2024 12 12 05;
    2024 12 12 06;
    2024 12 12 07;
    2024 12 12 08;
    2024 12 12 09;
    2024 12 12 10;
    2024 12 12 11;
    2024 12 12 12;
    2024 12 12 13;
    ];

index=index+1;

%% Redo Analysis
do_redo_analysis = 0;    % Do you want to run analysis on it?

if do_redo_analysis
    opts=struct;
    opts.do_ixon_main           = 0;   % ixon_main
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