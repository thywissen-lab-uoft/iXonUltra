clear composite_data
composite_data = struct;
index=1;

% Retake March data
% 200.9 G
composite_data(index).Name = '2024_12_13 200.9 G spectrum';
composite_data(index).Description = '2.5 Er, 200.9 G, 66.5 mW 4.5 Er pulse';
% composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [     
    2024 12 13 16;
    2024 12 13 17;
    2024 12 13 18;
    2024 12 13 19;
    2024 12 13 20;
    2024 12 13 21;
    2024 12 13 22;
    2024 12 14 01;
    2024 12 14 02;
    2024 12 14 03;
    ];
index = index +1;

% 190 G, 6.5 Er pulse
composite_data(index).Name = '2024_12_13 190 G spectrum';
composite_data(index).Description = '2.5 Er, 190 G, 66.5 mW 6.5 Er pulse';
% composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [     
    2024 12 14 04;
    2024 12 14 05;
    2024 12 14 06;
    2024 12 14 07;
    2024 12 14 08;
    2024 12 14 09;
    2024 12 14 10;
    2024 12 14 11;
    2024 12 14 12;
    2024 12 14 13;
    ];
index = index +1;

% 200.65 G, 5 Er pulse
composite_data(index).Name = '2024_12_13 200.65 G spectrum';
composite_data(index).Description = '2.5 Er, 200.65 G, 66.5 mW 5 Er pulse';
% composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [     
    2024 12 14 14;
    2024 12 14 15;
    2024 12 14 16;
    2024 12 14 17;
    2024 12 14 18;
    2024 12 14 19;
    2024 12 14 20;
    2024 12 14 21;
    2024 12 14 22;
    2024 12 14 23;
    

    ];
index = index +1;

% 195 G, 6 Er pulse
composite_data(index).Name = '2024_12_13 195 G spectrum';
composite_data(index).Description = '2.5 Er, 195 G, 66.5 mW 6 Er pulse';
% composite_data(index).Type = 'spectrum';
composite_data(index).Runs = [     
    2024 12 14 24;
    2024 12 14 25;
    2024 12 14 26;
    2024 12 14 27;
    2024 12 15 01;
    2024 12 15 02;
    2024 12 15 03;
    2024 12 15 04;
%     2024 12 15 05;no atom :-(
%     2024 12 15 06;

    ];
index = index +1;

% Overdirven

% % Retake retake March data
% % 200 G, 5 Er pulse
% composite_data(index).Name = '2024_12_19 200 G spectrum';
% composite_data(index).Description = '2.5 Er, 200 G, 66.5 mW 5 Er pulse';
% % composite_data(index).Type = 'spectrum';
% composite_data(index).Runs = [     
%     2024 12 19 01;
% %     2024 12 19 02; % no signal
%     2024 12 19 03;
%     2024 12 19 04;
%     2024 12 19 05;
%     2024 12 19 06;
%     2024 12 19 07;
%     2024 12 19 08;
% %     2024 12 20 01; % no signal
%     2024 12 20 02;
%     ];
% index = index +1;
% 
% % 200.9 G, 4.5 Er pulse
% composite_data(index).Name = '2024_12_19 200.9 G spectrum';
% composite_data(index).Description = '2.5 Er, 200.9 G, 66.5 mW 4.5 Er pulse';
% % composite_data(index).Type = 'spectrum';
% composite_data(index).Runs = [     
%     2024 12 20 03;
%     2024 12 20 04;
%     2024 12 20 05;
%     2024 12 20 06;
%     2024 12 20 07;
%     2024 12 20 08;
%     2024 12 20 09;
%     2024 12 20 10;
% %     2024 12 20 11; % no signal
%     2024 12 20 12;
%     ];
% index = index +1;

%%
composite_data = gatherCompositeData(composite_data);

%% Upload

doUpload = true;


GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
output_folder_name = '2024_12_13 Full spectrum vs U redo March';
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


