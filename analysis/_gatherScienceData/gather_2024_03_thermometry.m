clear composite_data
composite_data = struct;

%% Saving and Output

output_filename = '2024_03_thermometry';
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];


%% Define Runs Select runs
composite_data(1).Name = 'thermometry';
composite_data(1).Runs=[
    2024 03 01 12;
    2024 03 01 20;
    2024 03 04 04;
    2024 03 04 12;
    2024 03 05 05;
    2024 03 05 13;
    2024 03 06 06;
    2024 03 06 14;
    2024 03 07 5;
    2024 03 07 13;
    2024 03 08 04;
    2024 03 08 12;
    2024 03 14 05;
    2024 03 14 13;
    2024 03 19 05;
    2024 03 19 13;
    2024 03 25 05;
    2024 03 26 05;
    2024 03 26 13;
    2024 03 27 05;
    2024 03 27 13;
    2024 03 28 06;
    2024 03 28 14;
    2024 03 28 24;
    2024 03 28 32;
    ];


%% Gather All Data
opts=struct;
opts.MatFiles={'digdata.mat'};
composite_data = gatherCompositeData(composite_data,opts);

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