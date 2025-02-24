clear composite_data
composite_data = struct;

%%
runs = [  
    2025 01 27 03;
    2025 01 27 04;
    2025 01 27 05;
    2025 01 27 06;
    2025 01 27 07;
    2025 01 27 08;
    2025 01 27 09;
    2025 01 27 10;
    2025 01 28 01;
    2025 01 28 02;
    2025 01 28 03;
    2025 01 28 04;
    2025 01 28 05;
    2025 01 28 06;
    2025 01 28 07;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_01_27 200 G specta 65 mW, 6 Er pulse';
%%
runs = [  
    2025 01 28 08;
    2025 01 28 09;
    2025 01 28 10;
%     2025 01 28 11; %out of focus
    2025 01 28 12;
    2025 01 28 13;
    2025 01 28 14;
    2025 01 28 15;
    2025 01 28 16;
    2025 01 28 17;
    2025 01 28 18;
    2025 01 28 19;
    2025 01 28 20;
    2025 01 28 21;
    ];

composite_data(end+1).Runs = runs;
composite_data(end).Name = '2025_01_27 190 G specta 64.5 mW, 5.5 Er pulse';

composite_data(1)=[];
%% Gather All Data
composite_data = gatherCompositeData(composite_data);

%% Upload

doUpload = true;


GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name ='2025_01_27 Full Spectrum versus U';
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