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

    if ~isfield(opts,'doFit')
       opts.doFit = 1; 
    end

    if ~isfield(opts,'RadialStep')
        opts.RadialStep = 10;
    end
%% Compute Noise Threshold

noise_threshold = {};

for nn=1:length(bindata)
    noise_threshold{nn}=zeros(1,length(bindata(nn).LatticeBin));
    for jj=1:length(bindata(nn).LatticeBin)
            % site_area = bindata(nn).LatticeBin(1).lattice_spacing_px^2;
            site_area = bindata(nn).LatticeBin(1).LatticeSpacingPx^2;
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
            Zb0 = Zb;    

            % n1 = bindata(ii).LatticeBin(jj).n1;
            % n2 = bindata(ii).LatticeBin(jj).n2;
            n1 = 1:size(Zb,2);
            n2 = 1:size(Zb,1);

            % Compute Center
            [nn1, nn2]= meshgrid(n1, n2);
            binds = Zb<=noise_threshold{ii}(jj);
            Zb(binds)=0; % set sites below noise to zero
            Zb0(binds)=0;
            nn10 = nn1;
            nn20 = nn2;


            n2c=sum(nn2.*Zb,'all')/sum(Zb,'all');
            n1c=sum(nn1.*Zb,'all')/sum(Zb,'all');

            % Now remove all sites below noise
            allInds = 1:numel(Zb);
            allInds(binds)=[];
             Zb(binds)=[];
             nn2(binds)=[];
             nn1(binds)=[];

            [~,edges]=histcounts(Zb,opts.Bins);
            centers = (edges(2)-edges(1))+edges;
            centers(end)=[];    

            % Recenter the data
            X = nn1-n1c;
            Y = nn2-n2c;
            
            nn10 = nn10 -n1c;
            nn20 = nn20 - n2c;
            r0=abs(nn10+1i*nn20);




            % creating circular layers
            r_inds=round(abs(X+1i*Y)/opts.RadialStep)+1;
            
            r=accumarray(r_inds(:),abs(X(:)+1i*Y(:)),[],@mean);  
            n = accumarray(r_inds(:),Zb(:),[],@(x) numel(x));
            dev=accumarray(r_inds(:),Zb(:),[],@(x) std(x,1));

            data = zeros(numel(r),numel(centers));

            LatticeRadialHistogram(jj).RadialVector = r;
            LatticeRadialHistogram(jj).Edges = edges;
            LatticeRadialHistogram(jj).Centers = centers;
            LatticeRadialHistogram(jj).NPoints = n;
            LatticeRadialHistogram(jj).Sigma = dev/sqrt(n);

            for kk=1:max(r_inds)
                zR = Zb(r_inds(:)==kk);
                [n,edges]=histcounts(zR,edges);  
                data(kk,:)=n;

               
                
    
                cluster_number=2;
                try
                    [idx,c,sumD,D]=kmeans(zR(:),cluster_number);

                [c,inds]=sort(c,'ascend');
                sumD=sumD(inds);
                
                dd = zeros(numel(c),1);
                thresh = zeros(numel(c),1);

                % Zb0(allInds([r_inds(:)==kk]))=thresh(1);
                nCluster=zeros(numel(c),1);

                for bb=1:numel(c)
                    ind = inds(bb);
                    nThis = sum(idx==ind);
                    dd(bb) = sqrt(sumD(bb)/nThis);
                    
                    thresh(bb) = round(max(zR(idx==ind)));
                    nCluster(bb)=nThis;
                end
                thresh(end)=[];               
                
       
                out.ClusterNumber = cluster_number;
                out.ClusterThreshold = thresh;
                out.ClusterCentroid = c;
                out.ClusterRadius = dd;

                LatticeRadialHistogram(jj).ClusterNumber(kk) = nCluster(2);
                LatticeRadialHistogram(jj).ClusterCentroid(kk) = c(2);
                LatticeRadialHistogram(jj).ClusterThreshold(kk) = thresh(1);


                catch ME
                    LatticeRadialHistogram(jj).ClusterNumber(kk) = 0;
                    LatticeRadialHistogram(jj).ClusterCentroid(kk) = NaN;
                    LatticeRadialHistogram(jj).ClusterThreshold(kk) = 0;

                end
                
            end
            
            LatticeRadialHistogram(jj).N = data;  

            if opts.doFit
                myfit = fittype(@(A,s,r) A*exp(-r.^2/(2*s^2)),...
                    'coefficients',{'A','s'},'independent','r');
                fitopt = fitoptions(myfit);
                W=LatticeRadialHistogram(jj).ClusterNumber;
                W=W/max(W);
                W(1) = 0;
                fitopt.Weights=W;
                rThis = LatticeRadialHistogram(jj).RadialVector;
                rThis = rThis(:);
                tThis =  LatticeRadialHistogram(jj).ClusterThreshold;
                tThis = tThis(:);
                fitopt.StartPoint = [max(tThis) sum(tThis.*rThis)/sum(tThis)];
    
               fout = fit(rThis(:),tThis(:),myfit,fitopt);
    
                LatticeRadialHistogram(jj).GaussFit = fout;
                LatticeRadialHistogram(jj).Amplitude = fout.A;
                LatticeRadialHistogram(jj).Radius = fout.s;
                % LatticeRadialHistogram(jj).Background = fout.bg;

                thresh_map = feval(fout,r0);
                thresh_map = reshape(thresh_map,[size(Zb0,1) size(Zb0,2)]);
                thresh_min = 850;
                thresh_map(thresh_map<thresh_min)=thresh_min;
                LatticeRadialHistogram(jj).ZbinScale = Zb0./thresh_map;

                % 
                thresh_scale = feval(fout,r0)./fout.A;
                thresh_scale = reshape(thresh_scale,[size(Zb0,1) size(Zb0,2)]);
                thresh_min = 0.5;
                thresh_scale(thresh_scale<thresh_min)=thresh_min;
                LatticeRadialHistogram(jj).ZbinScale = Zb0./thresh_scale;

                % ZbinScaler
                
            end
       end
       bindata(ii).LatticeRadialHistogram=LatticeRadialHistogram;
   end    
end
