
%% Define Runs Select runs

runs=[
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
    clear digdata;


for kk=1:size(runs,1)
    dir_list = ixon_findRunDirectory(runs(kk,:));
    disp(repmat('-',1,60));disp(repmat('-',1,60));   
    disp(['loading ' num2str(runs(kk,:)) '...']);
  
    s1 = fullfile(imgdir,'Figures','digdata.mat');        
    %Load Files        
    data=load(s1);
    digdata(kk) = data;
end


%% Upload data
doUpload = 1;
GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24\shake_25Er'];

 try
     if ~exist(GDrive_root,'dir')
     mkdir(GDrive_root);
     end
 end
 
 if  doUpload && exist(GDrive_root,'dir')   
     disp('upload to google drive');
    gFile = fullfile(GDrive_root,'therm_data.mat');
    save(gFile,'digdata'); 
 end
