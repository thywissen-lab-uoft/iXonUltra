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
runs = [2023 12 14 06;
    2023 12 14 07;
    2023 12 14 08;
    2023 12 14 09;
    2023 12 14 10;
    2023 12 14 11;
    2023 12 14 12;
    2023 12 14 13;
    2023 12 14 14;
    2023 12 14 15;
    2023 12 14 16;
    2023 12 14 17;
    2023 12 14 18;
    2023 12 14 19;

    ];
%% Get the direcotry list
dir_list = ixon_findRunDirectory(runs);

for nn=1:length(dir_list)
   ixon_auto_dir = 0;
   imgdir = dir_list{nn};
   ixon_main;
    ixon_auto_dir = 1;
end



