clear composite_data
composite_data = struct;

%%
% 200.6 G
runs = [
    2025 03 15 15
    2025 03 15 16
    2025 03 15 17
    2025 03 15 18
    2025 03 15 19
    2025 03 15 20
    2025 03 15 21
    2025 03 15 22
    2025 03 15 23
    2025 03 15 24
    2025 03 15 25
    2025 03 15 26
    2025 03 15 27
    2025 03 15 28
    2025 03 15 29
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_15 200.6 G spectrum 54 mW, 2.5 Er pulse';

% 201.1 G
runs = [
    2025 03 15 31
    2025 03 15 32
    2025 03 16 01
    2025 03 16 02
    2025 03 16 03
    2025 03 16 04
    2025 03 16 05
    2025 03 16 06
    2025 03 16 07
    2025 03 16 08
    2025 03 16 09
    2025 03 16 10
    2025 03 16 11
    2025 03 16 12
    2025 03 16 13
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_16 201.1 G spectrum 54 mW, 3 Er pulse';

% 198.5 G
runs = [
    2025 03 16 15
    2025 03 16 16
    2025 03 16 17 
    2025 03 16 18
    2025 03 16 19
    2025 03 16 20
    2025 03 16 21
    2025 03 16 22
    2025 03 16 23
    2025 03 16 24
    2025 03 16 25
    2025 03 16 26
    2025 03 16 27
    2025 03 16 28
    2025 03 16 29
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_16 198.5 G spectrum 54.3 mW, 3 Er pulse';

% 200.35 G
runs = [
    2025 03 17 01
    2025 03 17 02
    2025 03 17 03 
    2025 03 17 04
    2025 03 17 05
    2025 03 17 06
    2025 03 17 07
    2025 03 17 08
    2025 03 17 09
    2025 03 17 10
    2025 03 17 11
    2025 03 17 12
    2025 03 17 13
    2025 03 17 14
    2025 03 17 15
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_17 200.35 G spectrum 54.3 mW, 3 Er pulse';

% 201.1 G 54.5mW, 3ER pulse
runs = [
    2025 03 17 17
    2025 03 17 18
    2025 03 17 19
    2025 03 17 20
    2025 03 17 21
    2025 03 17 22
    2025 03 17 23
    2025 03 17 24
    2025 03 17 25
    2025 03 17 26
    2025 03 17 27
    2025 03 17 28
    2025 03 17 29
    2025 03 17 30
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_17 201.1 G spectrum 54.5 mW, 3 Er pulse';

% 201.1 G 54 mW, 1 ER pulse
runs = [
    2025 03 18 04
    2025 03 18 05
    2025 03 18 06
    2025 03 18 07
    2025 03 18 08
    2025 03 18 09
    2025 03 18 10
    2025 03 18 11
    2025 03 18 12
    2025 03 18 13
    2025 03 18 14
    2025 03 18 15
    2025 03 18 16
    2025 03 18 17
    2025 03 18 18
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_18 201.1 G spectrum 54 mW, 1 Er pulse';

% 201.1 G 54.8 mW, 3 ER pulse
runs = [
    2025 03 18 20
    2025 03 18 21
    2025 03 18 22
    2025 03 18 23
    2025 03 18 24
    2025 03 18 25
    2025 03 18 26
    2025 03 18 27
    2025 03 18 28
    2025 03 18 29
    2025 03 18 30
    2025 03 19 01
    2025 03 19 02
    2025 03 19 03
    2025 03 19 04
    2025 03 19 01
    2025 03 19 02
    2025 03 19 03
    2025 03 19 04
    2025 03 19 06 %repeat
    2025 03 19 07 %repeat
    2025 03 19 08 %repeat
    2025 03 19 09 %repeat
    2025 03 19 10 %repeat
    2025 03 19 11 %repeat
    2025 03 19 12 %repeat
    ];
composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_03_18 201.1 G spectrum 54.8 mW, 3 Er pulse';

composite_data(1)=[];
%% Gather All Data
composite_data = gatherCompositeData(composite_data);

%% Upload

doUpload = true;


GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name ='2025_03_15 full spectrum';
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