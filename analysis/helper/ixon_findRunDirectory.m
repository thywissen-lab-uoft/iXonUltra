function [dirlist] = ixon_findRunDirectory(runs,src)

if nargin ==1
    src = 'X:\Data';
end

for kk=1:size(runs,1)
    % Construct strings for year, month, day, and run
    yStr = num2str(runs(kk,1));
    mStr = num2str(runs(kk,2),'%02d');
    dStr = num2str(runs(kk,3),'%02d');
    rStr = num2str(runs(kk,4),'%02d');

    % Find the location of the days data
    mDir = [yStr '.' mStr];
    dDir = [mStr '.' dStr];
    myDir = [yStr filesep mDir filesep dDir];
    myDirFull = fullfile(src,myDir);
    
    % Find all directories in this day
    myRuns = dir(myDirFull);    % Get folder contents    
    dirFlags = [myRuns.isdir];  % Flag the directories
    myRuns=myRuns(dirFlags);    % Get the directories
    myRuns = {myRuns.name};     % Get the names
    myRuns = myRuns(...         % Remove "fake" directories from dir
        ~ismember(myRuns ,{'.','..'}));
    
    % Get the first characteres to get the run number
    for nn=1:length(myRuns)
        s = myRuns{nn};
        if isequal(rStr,s(1:2))
           dirlist{kk} = fullfile(myDirFull,s); 
        end  
    end
    
end
end

