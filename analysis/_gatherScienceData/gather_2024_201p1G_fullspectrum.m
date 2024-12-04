clear composite_data
composite_data = struct;
index=1;
%% Introduction

%% 11/26-11/27
composite_data(index).Name = '2024_11_26 201.1 G 1.5 V Drive, Full Spectra, 70 mW Evap';
composite_data(index).Description = '2024_11_26 201.1 G 1.5 V Drive, Full Spectra, 70 mW Evap';
composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [     
    2024 11 26 12;
    2024 11 26 13;
    2024 11 26 14;
    2024 11 26 16;
    2024 11 26 16;
    2024 11 26 17;
    2024 11 26 18;
    2024 11 26 19;
    2024 11 27 01;
    2024 11 27 02;
    2024 11 27 03;
    2024 11 27 04;
    ];
index=index+1;

%% 11/27-11/28
composite_data(index).Name = '2024_11_27-28 201.1 G 2V and up Variable Drive, Full Spectra, 70 mW Evap';
composite_data(index).Description = '2024_11_27 201.1 G 2V and up Variable Drive, Full Spectra, 70 mW Evap';
composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [     
    2024 11 27 06;
    2024 11 27 07;
    2024 11 27 08;
    2024 11 27 09;
    2024 11 27 10;
    2024 11 28 01;
    2024 11 28 02;
    2024 11 28 03;
    ];
index=index+1;

%% 2024/11/28-2024/11/29
%201.1 G high field 50 ms mod ramp 4V drive, full spectrum
composite_data(index).Name = '11/28-11/29 201.1 G 4 V 70 mW';
composite_data(index).Description = '201.1 G high field, 2.5Er, 50 ms mod ramp, 4V, 70 mW evap,  Full Spectrum';
composite_data(index).Type = 'spectrum';
composite_data(index).Runs= [ 
        2024 11 28 04;
        2024 11 28 05;
        2024 11 28 06;
        2024 11 28 07;
        2024 11 28 08;
        2024 11 28 09;
        2024 11 28 10;
        2024 11 28 11;        
        2024 11 28 12;        
        2024 11 29 01;
        2024 11 29 02;
        2024 11 29 03;
        2024 11 29 04;
        2024 11 29 05;
        2024 11 29 06;
        2024 11 29 07;
        2024 11 29 08;
        2024 11 29 09;
        2024 11 29 10;
        2024 11 29 11;
        2024 11 29 12;
        2024 11 29 13;
        2024 11 29 14;
        2024 11 29 15;
        2024 11 29 16;
                
    ];
index=index+1;

%%

composite_data(index).Name = '12_01-02 201.1 G 4 V,full-spec,66.5 mW,5 Er pulse';
composite_data(index).Description = '12_01-02 201.1 G 4 V,full-spec,66.5 mW,5 Er pulse';
composite_data(index).Type = 'spectrum';
composite_data(index).Runs =[     
    2024 11 29 19;
    2024 11 29 20;
    2024 11 29 21;
    2024 11 29 22;
    2024 11 29 23;
    2024 11 29 24;
    2024 11 29 25;
    2024 11 29 26;
    2024 11 30 01;
    2024 11 30 02;
    2024 11 30 03;
    2024 11 30 04;
    2024 11 30 05;
    2024 11 30 06;
    2024 11 30 07;
    2024 11 30 08;
    2024 11 30 09;
    2024 11 30 10;
    2024 11 30 11;
    2024 11 30 12;
    2024 11 30 13;
    2024 11 30 14;
    2024 11 30 15;
    2024 11 30 16;
    2024 11 30 17;
    2024 11 30 18;

    ];
index=index+1;


%% 2024/12/01-2024/12/02
composite_data(index).Name = '2024_12_01-02 201.1 G 4 V Drive, Full Spectra, 67.5 mW Evap, 5 ER Pulse';
composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [     
    2024 12 01 02;
    2024 12 01 03;
    2024 12 01 04;
    2024 12 01 05;
    2024 12 01 06;
    2024 12 01 07;
    2024 12 01 08;
    2024 12 01 09;
    2024 12 01 10;
    2024 12 01 11;
    2024 12 01 12;
    2024 12 01 13;
    2024 12 01 14;
    2024 12 01 15;
    2024 12 01 16;
    2024 12 02 01;
    2024 12 02 02;
    2024 12 02 03;
    2024 12 02 04;
    2024 12 02 05;
    2024 12 02 06;
    2024 12 02 07;
    2024 12 02 08;
    2024 12 02 09;
    2024 12 02 10;
    2024 12 02 11;
    2024 12 02 12;
];
index=index+1;


%% Redo Analysis
do_redo_analysis = 0;    % Do you want to run analysis on it?

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


GDrive_root =['G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name = '2024_11_26 to 2024_12_02 Full Spectra High Temp';
saveDir = fullfile(GDrive_root,output_folder_name);

if doUpload
    try
        if ~exist(GDrive_root,'dir');mkdir(GDrive_root);end
        if ~exist(saveDir,'dir');mkdir(saveDir);end
         gFile = fullfile(saveDir,'composite_data.mat');

        disp(gFile);
        fprintf('upload to google drive ...');        
        save(gFile,'composite_data'); 
        disp('done!')   
    catch ME
        error('OH NO I COUD NOT UPL<OD');
    end
end 