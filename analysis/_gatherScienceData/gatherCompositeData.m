function composite_data=gatherCompositeData(composite_data,opts)
% gatherBulk_digdata.m
%
% composite_data : A structure array
% composite_data(index).Name;
% composite_data(index.Runs;
%
% opts : option structure;
% opts.MatFiles : which analysis to do

if nargin ==1
    opts=struct;
end

if ~isfield(opts,'MatFiles')
    opts.MatFiles = {'conductivity_data.mat','digdata.mat'};
end

for kk=1:length(composite_data)    
    dir_list = ixon_findRunDirectory(composite_data(kk).Runs);
    disp(repmat('-',1,60));
    disp(['Bulk load ' composite_data(kk).Name ' ' num2str(length(dir_list)) ' runs.']);
    for nn=1:length(dir_list)
        imgdir = dir_list{nn};
        for jj=1:length(opts.MatFiles)
            filename = fullfile(imgdir,'Figures',opts.MatFiles{jj});
            [~,field_name,~]=fileparts(filename);
            composite_data(kk).(field_name)(nn) = load(filename);
        end
    end


end
