function saveCompositeData(composite_data,opts)

if nargin ==1
   opts=struct; 
end

if ~isfield(opts,'doUpload')
    opts.doUpload = false;    
end

if ~isfield(opts,'Name')
   opts.Name = 'composite_analysis'; 
end

%%

GDrive_root =['G:\My Drive\Lattice Shared\SharedData\Conductivity_Saturated_23-24'];
Local_root ='X:\SaturatedConductivity';

%%
try
    if ~exist(Local_root,'dir')
        mkdir(Local_root); 
    end   
    local_target = fullfile(Local_root,opts.Name);
    
    if ~exist(local_target,'dir')
        mkdir(local_target);
    end
    
    fprintf('saving locally ...');
    myfile = fullfile(local_target,'composite_data.mat');
    save(myfile,'composite_data'); 
    disp('done!')
catch ME
    keyboard
end


%% Upload the data
if opts.doUpload    
     try
         if ~exist(GDrive_root,'dir')
         mkdir(GDrive_root);
         end
         
         gdrive_target = fullfile(GDrive_root,opts.Name);
         if ~exist(gdrive_target,'dir')
         mkdir(gdrive_target);
         end

        fprintf('uploading to google drive ...');
        myfile = fullfile(gdrive_target,'composite_data.mat');
        save(myfile,'composite_data'); 
        disp('done!')
         
     catch ME
         keyboard
     end    
end

