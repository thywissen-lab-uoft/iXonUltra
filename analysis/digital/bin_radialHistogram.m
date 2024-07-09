function bindata = bin_radialHistogram(bindata,opts)
    if ~isfield(bindata,'LatticeBin')
       return;
    end   
    
%% Options
    if nargin == 1
        opts = struct;
    end
    
    if ~isfield(opts,'ROI')
        opts.ROI = 'max';
    end   
   
    if ~isfield(opts,'Nthresh')
       opts.Nthresh = 1000; 
    end
    
    if ~isfield(opts,'Bins')
        opts.Bins = 30;
    end   
    
    if ~isfield(opts,'doAnimate')
       opts.doAnimate = 0; 
    end

    if ~isfield(opts,'RadialStep')
        opts.RadialStep = 10;
    end
%% Compute Noise Threshold

noise_threshold = {};

for nn=1:length(bindata)
    noise_threshold{nn}=zeros(1,length(bindata(nn).LatticeBin));
    for jj=1:length(bindata(nn).LatticeBin)
            site_area = bindata(nn).LatticeBin(1).lattice_spacing_px^2;
            noise_per_px=bindata(nn).NoiseEstimation(1);
            noise_variance_per_px = noise_per_px.^2;
            noise_variance_per_site = site_area*noise_variance_per_px;
            noise_per_site = sqrt(noise_variance_per_site);
            noise_threshold_this = 2*noise_per_site;
            noise_threshold{nn}(jj)=max([noise_threshold_this 0]);
    end
end
   
 %% Compute Center Points
   for ii = 1:length(bindata)
       LatticeRadialHistogram = struct;
       for jj=1:length(bindata.LatticeBin)    
            % Get Data
            Zb = bindata(ii).LatticeBin(jj).Zbin;    
            % n1 = bindata(ii).LatticeBin(jj).n1;
            % n2 = bindata(ii).LatticeBin(jj).n2;
            n1 = 1:size(Zb,2);
            n2 = 1:size(Zb,1);

            % Compute Center
            [nn1, nn2]= meshgrid(n1, n2);
            binds = Zb<=noise_threshold{ii}(jj);
            Zb(binds)=0; % set sites below noise to zero

            n2c=sum(nn2.*Zb,'all')/sum(Zb,'all');
            n1c=sum(nn1.*Zb,'all')/sum(Zb,'all');

            % Now remove all sites below noise
            Zb(binds)=[];
            nn2(binds)=[];
            nn1(binds)=[];

            [~,edges]=histcounts(Zb,opts.Bins);
            centers = (edges(2)-edges(1))+edges;
            centers(end)=[];    

            % Recenter the data
            X = nn1-n1c;
            Y = nn2-n2c;

            % creating circular layers
            r_inds=round(abs(X+1i*Y)/opts.RadialStep)+1;
            r=accumarray(r_inds(:),abs(X(:)+1i*Y(:)),[],@mean);  
            n = accumarray(r_inds(:),Zb(:),[],@(x) numel(x));

            data = zeros(numel(r),numel(centers));

            LatticeRadialHistogram(jj).RadialVector = r;
            LatticeRadialHistogram(jj).Edges = edges;
            LatticeRadialHistogram(jj).Center = centers;
            LatticeRadialHistogram(jj).NPoints = n;

            for kk=1:max(r_inds)
                zR = Zb(r_inds(:)==kk);
                [n,edges]=histcounts(zR,edges);  
                data(kk,:)=n;
            end
            LatticeRadialHistogram(jj).N = data;  
       end
       bindata(ii).LatticeRadialHistogram=LatticeRadialHistogram;
   end    
end
