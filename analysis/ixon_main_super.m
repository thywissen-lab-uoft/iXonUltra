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

runs = [2024 01 16 23;
    2024 01 16 24;
    2024 01 16 25;
    2024 01 16 26;
    2024 01 16 27;
    2024 01 16 28;
    2024 01 16 29;
    2024 01 16 30;
    ];

%% Get the direcotry list
dir_list = ixon_findRunDirectory(runs);

for nn=1:length(dir_list)
   ixon_auto_dir = 0;
   imgdir = dir_list{nn};
   ixon_main;
    ixon_auto_dir = 1;
end



