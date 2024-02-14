% ixon_main_super.m
%
% Author : CF Fujiwara
%

% Display this filename
disp(repmat('-',1,60));disp(repmat('-',1,60));    
disp(['Calling ' mfilename '.m']);
disp(repmat('-',1,60));disp(repmat('-',1,60));    

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath));

a = fileparts(curpath);
addpath(a);addpath(genpath(a));

src = 'X:\Data';

%%
% runs = [2023 12 14 06;
%     2023 12 14 07;
%     2023 12 14 08;
%     2023 12 14 09;
%     2023 12 14 10;
%     2023 12 14 11;
%     2023 12 14 12;
%     2023 12 14 13;
%     2023 12 14 14;
%     2023 12 14 15;
%     2023 12 14 16;
%     2023 12 14 17;
%     2023 12 14 18;
%     2023 12 14 19;
% 
%     ];

% runs = [2024 01 03 03;
%     2024 01 03 04;
%     2024 01 03 05;
%     2024 01 03 06;
%     2024 01 03 07;
%     2024 01 03 08;
%     2024 01 03 09;
%     2024 01 03 10;
%     2024 01 03 11;
%     2024 01 03 12;
%     2024 01 03 13;
%     ];

runs = [2024 02 13 20;
2024 02 13 21;
2024 02 13 22;    
2024 02 13 23;    
2024 02 13 24;    
2024 02 13 25;    
2024 02 13 26;    
2024 02 13 27;    
2024 02 13 28;    
2024 02 13 29;    
2024 02 13 30;    
2024 02 13 31;    
    ];

%% Get the direcotry list
dir_list = ixon_findRunDirectory(runs);

for nn=1:length(dir_list)
   ixon_auto_dir = 0;
   imgdir = dir_list{nn};
   ixon_main;
    ixon_auto_dir = 1;
end

%% Get the direcotry list
dir_list = ixon_findRunDirectory(runs);

for nn=1:length(dir_list)
    bin_auto_file = 0;
    filename = fullfile(dir_list{nn},'figures','bindata.mat');
    bin_imgdir = fullfile(dir_list{nn},'figures');
    ixon_bin_analysis;
    bin_auto_file = 1;    
end



%% Get the direcotry list
dir_list = ixon_findRunDirectory(runs);

for nn=1:length(dir_list)
    dig_auto_file = 0;
    dig_imgdir = fullfile(dir_list{nn},'figures');
    filename = fullfile(dir_list{nn},'figures','digdata.mat');

    ixon_dig_analysis;
    dig_auto_file = 1;    
end

