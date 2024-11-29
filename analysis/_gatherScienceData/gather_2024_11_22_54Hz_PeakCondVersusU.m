clear composite_data
composite_data = struct;
index=1;


%% 2024/11/22-2024/11/23
%201.1 G high field 50 ms mod ramp 11/22-11/23 vary evap depth
composite_data(index).Name = '11/22-11/23 201.1 G 54 Hz';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 54 Hz 0.4V, vary evap depth ';
composite_data(index).Runs= [ 
        2024 11 22 23;
        2024 11 22 24;
        2024 11 22 25;
        2024 11 22 26;
        2024 11 23 01;
        2024 11 23 03;
        2024 11 23 04;    
        2024 11 23 09;
        2024 11 23 10;
        2024 11 23 11;
        2024 11 23 12;
        2024 11 23 14;
        2024 11 23 16;
        2024 11 23 17;
        2024 11 23 19; % questionable
    ];
index=index+1;


%% 2024/11/23-2024/11/24
%201.1 G high field 50 ms mod ramp 11/23-11/24
composite_data(index).Name = '11/23-11/24 201.1 G, 54 Hz, vary evap and 4 Er pulse';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 54 Hz 0.4V, vary evap depth, pulse 4 ER lattice ';
composite_data(index).Runs= [ 
        2024 11 23 20; % atom number fluct affecting fit, questionalbe SNR
        2024 11 23 21; % atom number drifting 
%         2024 11 23 22; % bad SNR
%         2024 11 23 23; % bad SNR
%         2024 11 23 24; % bad SNR
%         2024 11 23 25; % bad SNR
%         2024 11 23 26; % bad SNR
        2024 11 23 27;
        2024 11 23 28;
        2024 11 24 01;
%         2024 11 24 02; % questionable fit
%         2024 11 24 03; % bad SNR
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

%% 2024/11/25-2024/11/26
%201.1 G high field 50 ms mod ramp 2024/11/25-2024/11/26 vary evap depth
composite_data(index).Name = '54 Hz 201.1 G 0.8V vary evap depth';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 54 Hz 0.8V, vary evap depth ';
composite_data(index).Runs =[     
    2024 11 25 22;
    2024 11 25 23;% questionable
    2024 11 25 24; % bad SNR
    2024 11 26 01;
    2024 11 26 02;
    2024 11 26 03;% good 68.5 mW
    2024 11 26 04;% atom number varies
    2024 11 26 05;
    2024 11 26 06;% questionable fit - below 0.5um, 69.5 mW
    2024 11 26 07;% questionable fit
    2024 11 26 08;% questionable fit 
    2024 11 26 09; % questionable fit 70 mW
    2024 11 26 10;
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
%% Upload

doUpload = true;


GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name = '2024_11_20 Peak Cond 54 Hz Versus U';
saveDir = fullfile(GDrive_root,output_folder_name);

if doUpload
    try
        if ~exist(GDrive_root,'dir');mkdir(GDrive_root);end
        if ~exist(GDrive_root,'dir');mkdir(saveDir);end
         gFile = fullfile(saveDir,'composite_data.mat');

        disp(gFile);
        fprintf('upload to google drive ...');        
        save(gFile,'composite_data'); 
        disp('done!')   
    catch ME
        error('OH NO I COUD NOT UPL<OD');
    end
end 

