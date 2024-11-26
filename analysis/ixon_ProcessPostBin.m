function data = ixon_ProcessPostBin(data,opts)
%% Options
if nargin == 1
    opts = struct;
    
    opts.CompensateMethod           ='gauss';'none';'custom';
    % If radially compensating
    opts.CompensateGaussRadius      = 50;
    opts.CompensateMax              = 1.5;
    % If using a custom map
    opts.CompensateCustomMap        = 'asdfassdf.mat';
    opts.ClusterNumber              =2;
end

%% Error Checking
if ~isfield(data,'LatticeBin')
    error('no LatticeBin')
end

%% Recenter bindata
if length(data)>1
    data = bin_recenter(data);
end

%% Preparing
for kk=1:length(data)
    for rr = 1:length(data(kk).LatticeBin)
        fprintf(['(' num2str(kk) '/' num2str(length(data)) ')' num2str(rr) ' '])

        %% Spatial Compensation
        fprintf('spatial compensate...');
        if ~isfield(data(kk).LatticeBin(rr),'ZbinRaw')
            ZbinRaw = data(kk).LatticeBin(rr).Zbin;
        else
            ZbinRaw = data(kk).LatticeBin(rr).ZbinRaw;
        end

        switch opts.CompensateMethod
            case 'none'
                Zbin = ZbinRaw;
            case 'gauss'
                n1 = data(kk).LatticeBin(rr).n1;
                n2 = data(kk).LatticeBin(rr).n2;
                [nn1,nn2]=meshgrid(n1,n2);
                Ztemp = ZbinRaw;
                Ztemp(isnan(Ztemp)) = 0;
                Ztemp(isinf(Ztemp)) = 0;
                
                n1c=sum(sum(Ztemp,3).*nn1,'all')/sum(Ztemp,'all');
                n2c=sum(sum(Ztemp,3).*nn2,'all')/sum(Ztemp,'all');
                s=opts.CompensateGaussRadius;   
                N = opts.CompensateMax;
                map = exp(-(nn1-n1c).^2/(2*s^2)).*exp(-(nn2-n2c).^2/(2*s^2));
                map_inv = 1./map;
                map_inv(map_inv>N)=N;
                Zbin = ZbinRaw.*map_inv;   
            case 'custom'
                % Load a file which has a pixel map
            case 'otherwise'
                error('invalid spatial CompensateMethod')
        end        
        %% Clustering
        fprintf('kmeans...');

        z=Zbin(:);z(z==0)=[];       % Get histogram of all counts
        z(isnan(z))=[];
        % Sort into clusters
        [idx,ClusterCentroids,sumD,D]=kmeans(z,opts.ClusterNumber);
        
        % Sort by centroid
        [ClusterCentroids,inds]=sort(ClusterCentroids,'ascend');
        sumD=sumD(inds);
        
        ClusterRadii = zeros(numel(ClusterCentroids),1);
        ClusterThresholds = zeros(numel(ClusterCentroids),1);
        for jj=1:numel(ClusterCentroids)
            ind = inds(jj);
            nThis = sum(idx==ind);
            ClusterRadii(jj) = sqrt(sumD(jj)/nThis);
            ClusterThresholds(jj) = round(max(z(idx==ind)));
        end
        ClusterThresholds(end)=[];    

        %% Fit PDF To Cluster
        fprintf('pdf...');

        boundLow = ClusterCentroids(2)-1.5*ClusterRadii(2);% send this to be based on cluster threshold
        boundLow = ClusterThresholds(1);
        boundHigh = 2*ClusterCentroids(2);
        boundHigh = inf;
        
        z_truncate =z;
        z_truncate(z_truncate<boundLow)=[];
        z_truncate(z_truncate>boundHigh)=[];

        [pdf1_c,pdf1_cint] = mle(z_truncate,'distribution','normal',...
            'TruncationBounds',[boundLow boundHigh]);
        
     
        
        PDF1_Distribution = 'normal';
        PDF1_Center = pdf1_c(1);
        PDF1_Radius = pdf1_c(2);
        PDF1_Center_Bounds = pdf1_cint(:,1);
        PDF1_Radius_Bounds = pdf1_cint(:,2);
        PDF1 = @(x) pdf('normal',x,PDF1_Center,PDF1_Radius);

        data(kk).LatticeBin(rr).PostBinOptions=opts;
        data(kk).LatticeBin(rr).Zbin=Zbin;
        data(kk).LatticeBin(rr).ClusterThreshold     = ClusterThresholds;
        data(kk).LatticeBin(rr).ClusterCentroids     = ClusterCentroids;
        data(kk).LatticeBin(rr).ClusterRadii         = ClusterRadii;
        data(kk).LatticeBin(rr).PDF1                 = PDF1;
        data(kk).LatticeBin(rr).PDF1_Distribution    = PDF1_Distribution;
        data(kk).LatticeBin(rr).PDF1_Center          = PDF1_Center;
        data(kk).LatticeBin(rr).PDF1_Radius          = PDF1_Radius;
        data(kk).LatticeBin(rr).PDF1_Center_Bounds   = PDF1_Center_Bounds;
        data(kk).LatticeBin(rr).PDF1_Radius_Bounds   = PDF1_Radius_Bounds;
        disp('done');
    end
end
end

