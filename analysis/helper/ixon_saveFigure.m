function ixon_saveFigure(atomdata,hF,filename)
global ixon_imgdir
ext='.png';
save_qual='-r90';


% The directory where the figures are saved
figDir=fullfile(ixon_imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

% Create the name of the figure
[filepath,name,~]=fileparts(ixon_imgdir);


% Make the figure name with the location
saveLocation=fullfile(figDir,[filename ext]);
% saveLocation='C:

% Save the figure and the png
fprintf([datestr(now,13) ' Saving figure handle to ']);
fprintf([filename ext ' ... ']);
set(0,'CurrentFigure', hF);
set(hF,'PaperPositionMode','auto');
print('-dpng',save_qual,saveLocation);
disp('Saved!');

% savefig(the_figure,saveLocation);
end

