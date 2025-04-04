function out = bin_NormalizedStatistics(bindata)

out=struct;
%% Init Data
Npics = length(bindata(1).LatticeBin);
count_center = zeros(length(bindata),Npics);
count_sigma  = zeros(length(bindata),Npics);

%% Get all Count Centers
for kk=1:length(bindata)
    for rr=1:length(bindata(kk).LatticeBin)
        if isfield(bindata(kk).LatticeBin(rr),'PDF1_Center') && ~isempty(bindata(kk).LatticeBin(rr).PDF1_Center)
       
            count_center(kk,rr) = bindata(kk).LatticeBin(rr).PDF1_Center;
            count_sigma(kk,rr) = bindata(kk).LatticeBin(rr).PDF1_Radius; 
        else
            count_center(kk,rr) = NaN;
            count_sigma(kk,rr) = NaN;
        end
    end
end

    %% Get all Z bin data
raw={};
normalized={};
processed={};
for rr=1:length(bindata(1).LatticeBin)
    ZRaw = zeros(size(bindata(1).LatticeBin(1).Zbin,1),...
        size(bindata(1).LatticeBin(1).Zbin,2),length(bindata));
    ZProcessed = ZRaw;
    ZNormalized =ZRaw;
    
    for kk=1:length(bindata)
        ZRaw(:,:,kk)= bindata(kk).LatticeBin(rr).ZbinRaw;
        ZProcessed(:,:,kk) = bindata(kk).LatticeBin(rr).Zbin;   
        center_dev = (count_center(kk,rr)-median(count_center(:,rr)))/median(count_center(:,rr));
        if abs(center_dev) > 0.75
            C = median(count_center(:,rr));
            count_center(kk,rr)=C;
            count_sigma(kk,rr) = median(count_sigma(:,rr));
        else
            C = count_center(kk,rr);
        end        
        ZNormalized(:,:,kk) = bindata(kk).LatticeBin(rr).Zbin/C;
    end
    
    raw{rr}=ZRaw;
    normalized{rr}=ZNormalized;
    processed{rr}=ZProcessed;
end

%% 

for rr=1:length(bindata(1).LatticeBin)
    z = normalized{rr};
    zall=z(:);
    zall(zall==0)=[];
    zall(isnan(zall))=[];

    zall_high=zall;
    zall_high(zall<0.5) = [];
    [f,xf] = ksdensity(zall_high);
    [~,ind]=max(f);
    n0=xf(ind);

    zall_high(zall_high<n0)=[];

    
    pd=fitdist(zall_high-n0,'Half Normal');

    pd = makedist('Normal',n0,pd.sigma);
    out(rr).PDF1                = pd;

    out(rr).PDF1_Distribution   = 'Normal';
    out(rr).PDF1_Center         = pd.mu;
    out(rr).PDF1_Radius         = pd.sigma;
    
    % 
    %    % Sort into clusters
    % [idx,ClusterCentroids,sumD,D]=kmeans(zall,2);
    % 
    % % Sort by centroid
    % [ClusterCentroids,inds]=sort(ClusterCentroids,'ascend');
    % sumD=sumD(inds);
    % 
    % ClusterRadii = zeros(numel(ClusterCentroids),1);
    % ClusterThresholds = zeros(numel(ClusterCentroids),1);
    % for jj=1:numel(ClusterCentroids)
    %     ind = inds(jj);
    %     nThis = sum(idx==ind);
    %     ClusterRadii(jj) = sqrt(sumD(jj)/nThis);
    %     ClusterThresholds(jj) = max(zall(idx==ind));
    % end
    % ClusterThresholds(end)=[];     
    % 
    % boundLow  = ClusterThresholds(1);
    % boundHigh = 2*ClusterCentroids(2);
    % 
    % 
    % boundLow  = 0.9*ClusterCentroids(2);
    % boundHigh = 3.0*ClusterCentroids(2);
    % 
    % 
    % z_truncate =zall;
    % z_truncate(z_truncate<boundLow)=[];
    % z_truncate(z_truncate>boundHigh)=[];
    % 
    % [pdf1_c,pdf1_cint] = mle(z_truncate,'distribution','normal',...
    %     'TruncationBounds',[boundLow boundHigh]);
    % 
    % PDF1_Distribution = 'normal';
    % PDF1_Center = pdf1_c(1);
    % PDF1_Radius = pdf1_c(2);
    % PDF1_Center_Bounds = pdf1_cint(:,1);
    % PDF1_Radius_Bounds = pdf1_cint(:,2);
    % PDF1 = @(x) pdf('normal',x,PDF1_Center,PDF1_Radius);
        
    % out(rr).ClusterThreshold    = ClusterThresholds;
    % out(rr).ClusterCentroids    = ClusterCentroids;
    % out(rr).ClusterRadii        = ClusterRadii;    
    % out(rr).PDF1                = PDF1;
    % out(rr).PDF1_Distribution   = PDF1_Distribution;
    % out(rr).PDF1_Center         = PDF1_Center;
    % out(rr).PDF1_Radius         = PDF1_Radius;
    % out(rr).PDF1_Center_Bounds  = PDF1_Center_Bounds;
    % out(rr).PDF1_Radius_Bounds  = PDF1_Radius_Bounds;

end

end

