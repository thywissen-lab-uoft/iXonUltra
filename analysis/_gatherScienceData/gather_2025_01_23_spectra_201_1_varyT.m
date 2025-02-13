clear composite_data
composite_data = struct;

%%
% 64 mW, 4 Er pulse
composite_data(end+1).Name = '2025_01_23 201.1G specta 64 mW 4 Er pulse';
composite_data(end).Runs = [
    2025 01 23 07;
    2025 01 23 08;
    2025 01 23 09;
    2025 01 23 10;
    2025 01 23 11;
    2025 01 23 12;
    2025 01 24 01;
    2025 01 24 02;
    2025 01 24 03;
    2025 01 24 04;
    2025 01 24 05;
    2025 01 24 06;
    2025 01 24 07;
    2025 01 24 08;
    2025 01 24 09;
    2025 01 24 10; %repeat
    2025 01 24 11; %repeat
    2025 01 24 12; %repeat
    ];

% 66 mW, 5Er pulse
composite_data(end+1).Name = '2025_01_24 201.1G specta 66 mW 5 Er pulse';
composite_data(end).Runs = [
    2025 01 24 15;
    2025 01 24 16;
    2025 01 24 17;
    2025 01 24 18;
%     2025 01 24 19; % out of focus
    2025 01 24 20;
    2025 01 24 21;
    2025 01 24 22;
    2025 01 25 01;
    2025 01 25 02;
    2025 01 25 03;
    2025 01 25 04;
    2025 01 25 05;
    2025 01 25 06;
    2025 01 25 07;
%     2025 01 25 08; % out of focus
    2025 01 25 09;
    2025 01 25 10;
    2025 01 25 11;
    2025 01 25 12;
    2025 01 25 13;
    2025 01 25 14;
    2025 01 25 15;
    2025 01 25 16;
    2025 01 25 17;
    2025 01 25 18;
    2025 01 25 19;
    2025 01 25 20;
    2025 01 25 21;
    2025 01 25 22;
%     2025 01 25 23; out of focus
    ];

% 63.8 mW
composite_data(end+1).Name = '2025_01_25 201.1G specta 63.8 mW';
composite_data(end).Runs = [
    2025 01 25 25;
    2025 01 25 26;
    2025 01 25 27;
    2025 01 26 01;
    2025 01 26 02;
    2025 01 26 03;
    2025 01 26 04;
    2025 01 26 05;
    2025 01 26 06;
    2025 01 26 07;
    2025 01 26 08;
    2025 01 26 09;
    2025 01 26 10;
    2025 01 26 11;
    2025 01 26 12;
    2025 01 26 13;
    2025 01 26 14;
    2025 01 26 15;
    2025 01 26 16;
    2025 01 26 17;
    2025 01 26 18;
    2025 01 26 19;
    2025 01 26 20;
%     2025 01 26 21; % out of focus
%     2025 01 26 22; % out of focus
%     2025 01 26 23; % out of focus
    2025 01 26 24;
    2025 01 26 25;
    2025 01 26 26;
    ];

composite_data(1)=[];
%% Gather All Data
composite_data = gatherCompositeData(composite_data);

%% Upload

doUpload = true;


GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name ='2025_01_23 201_1 G Full Spectrum versus T';
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