clear composite_data
composite_data = struct;
index=1;
%% 2024/12/02-2024/12/03

%Single frequency 54 Hz, 0.4 V, Variable Field
% 2.5 Er, 50 ms mod ramp, 68 mW
composite_data(index).Name = '12/02 0.4V 68 mW 54 Hz';
composite_data(index).Description = '12/02 0.4V 68 mW 54 Hz';
composite_data(index).Type = 'peak';
composite_data(index).Runs= [ 
        2024 12 02 16;
        2024 12 02 18
        2024 12 02 20;
        2024 12 03 01;
        2024 12 03 03;
        2024 12 03 05;
        2024 12 03 07;
        2024 12 03 09;
        2024 12 03 11;
        2024 12 03 13;
        2024 12 03 15;

    ];
index=index+1;
%% 2024/12/02-2024/12/03

%Single frequency 54 Hz, 0.8 V, Variable Field
% 2.5 Er, 50 ms mod ramp, 68 mW
composite_data(index).Name = '12/02 0.8V 68 mW 54 Hz';
composite_data(index).Description = '12/02 0.8V 68 mW 54 Hz';
composite_data(index).Type = 'peak';
composite_data(index).Runs= [ 
        2024 12 02 17;
        2024 12 02 19
        2024 12 02 21;
        2024 12 03 02;
        2024 12 03 04;
        2024 12 03 06;
        2024 12 03 08;
        2024 12 03 10;
        2024 12 03 12;
        2024 12 03 14;
        2024 12 03 16;
    ];
index=index+1;

%% Redo Analysis
do_redo_analysis =0;    % Do you want to run analysis on it?

if do_redo_analysis
    opts=struct;
    opts.do_ixon_main           = 0;   % ixon_main
    opts.do_ixon_bin_analysis   = 0;   % ixon_bing
    opts.do_ixon_dig_analysis   = 1;   % ixon_dig
    ixon_super(composite_data,opts)
end

%% Gather Data
composite_data = gatherCompositeData(composite_data);

%% Upload

doUpload = true;


GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name ='2024_12_02 54Hz peak Cond variable field';
saveDir = fullfile(GDrive_root,output_folder_name);

if doUpload
    try
        if ~exist(GDrive_root,'dir');mkdir(GDrive_root);end
        if ~exist(saveDir,'dir');mkdir(saveDir); end
         gFile = fullfile(saveDir,'composite_data.mat');

        disp(gFile);
        fprintf('upload to google drive ...');        
        save(gFile,'composite_data'); 
        disp('done!')   
    catch ME
        error('OH NO I COUD NOT UPL<OD');
    end
end 

