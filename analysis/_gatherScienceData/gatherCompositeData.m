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
    opts.MatFiles = {'conductivity_data.mat','digdata.mat','dig_radial_data.mat','qpd.mat'};
%         opts.MatFiles = {'conductivity_data.mat','digdata.mat','dig_radial_data.mat'};

    %opts.MatFiles = {'conductivity_data.mat','digdata.mat'};
end

% for kk=1:length(opts.MatFiles)
%     [~,field_name,~]=fileparts(opts.MatFiles{kk});
%     [composite_data.(field_name)]=deal(struct);
% end

for kk=1:length(composite_data)    
    dir_list = ixon_findRunDirectory(composite_data(kk).Runs);
    disp(repmat('-',1,60));
    disp(['Bulk load ' composite_data(kk).Name ' ' num2str(length(dir_list)) ' runs.']);
    for nn=1:length(dir_list)
        imgdir = dir_list{nn};
%         keyboard
        for jj=1:length(opts.MatFiles)
            filename = fullfile(imgdir,'Figures',opts.MatFiles{jj});
            [~,field_name,~]=fileparts(filename);
%             keyboard
            if exist(filename)
                composite_data(kk).(field_name)(nn) = load(filename);
            else
                warning(imgdir);
                warning('no file found')
                composite_data(kk).(field_name)(nn) = [];
            end
        end
    end


end

