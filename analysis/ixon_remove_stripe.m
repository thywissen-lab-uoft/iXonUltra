function ixon_remove_stripe


% Choose the directory where the images to analyze are stored
disp([datestr(now,13) ' Choose an image analysis folder...']);
dialog_title='Choose the root directory of the images';
default_analysis_dir = ixon_getDayDir;


newdir=uigetdir(default_analysis_dir,dialog_title);

names=dir([newdir filesep 'iXon*.mat']);
names={names.name};
    
[a,b,c]=fileparts(newdir);

outdir = [a filesep 'stripes ' b c];

mkdir(outdir);
    
for mm=1:length(names)
   filename=fullfile(a,[b c],names{mm});
   d=load(filename);
   data=d.data;
   F = [data.Flags];
   
   if F.plane_selection_dotilt
       disp(['Moving' num2str(filename)]);
       movefile(filename,outdir);
   end
end

end

