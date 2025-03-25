function  [figs,out] = bin_compareHistograms(bindata,opts)
% This function accounts for 


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

% Remove NaN Data (ie. bad counts/binning)
count_center(isnan(count_center))=[];
count_sigma(isnan(count_sigma))=[];

% Calculate histogram max and dividing line for plotting
histMax = median(count_center(:)) + 3*median(count_sigma(:));
histMax = ceil(histMax/1e3)*1e3;
histDivider = median(count_center(:))- 2*median(count_sigma(:));

% Compute Histgoram Edges, Center, and Indeces
edges=linspace(0,histMax,50);
centers = (edges(1:end-1) + edges(2:end))/2;    
iL = centers<=histDivider;
iH = ~iL; 

% Calculate the same for the normalized image
edges_norm = edges/mean(count_center(:));
histDividerNorm = histDivider/mean(count_center(:));

centers_norm = (edges_norm(1:end-1) + edges_norm(2:end))/2;    
iL_norm = centers_norm<=histDividerNorm;
iH_norm = ~iL_norm; 



%% Normalized Fit

out = bin_NormalizedStatistics(bindata);
%%
ca = [0 0 0];       
cb = [0.7 .1 .6];
cc = [linspace(ca(1),cb(1),1000)' ...
    linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
        
for rr=1:length(bindata(1).LatticeBin)
    n1 = bindata(1).LatticeBin(rr).n1;
    n2 = bindata(1).LatticeBin(rr).n2;

    
    figs(rr) = figure;
    figs(rr).Color='w';
    figs(rr).Position=[50 50 1600 300];  
    colormap(figs(rr),cc);        
    hold on
    
    % Compute Raw Histogram
    zr=raw{rr}(:);
    zr(zr==0)=[];
    zr(isnan(zr))=[];    
    rawN = histcounts(zr,edges); 

    % Compute Spatial Compensate Histogram
    zp = processed{rr}(:);
    zp(zp==0)=[];
    zp(isnan(zp))=[];
    spatialN = histcounts(zp,edges);

    % Compute Spatial Compensate + Normalized Histogram
    zn=normalized{rr}(:);
    zn(zn==0)=[];
    normN = histcounts(zn,edges_norm); 


    % Raw Histgoram
    subplot(141);
    pHistB1 = bar(centers(iL),rawN(iL),'linestyle','none',...
        'facecolor','k');
    xlim([0 max(edges)]);    
    ylabel('occurences');
    xlabel('counts per lattice site');
    ylim([0 rawN(2)]);
    hold on
    % High Counts
    yyaxis right
    pHistB2 = bar(centers(iH),rawN(iH),'linestyle','none',...
    'FaceColor',[0.6 0 0.5]);
    set(gca,'box','on','YColor',[0.6 0 0.5]*.8);
    ylabel('occurences');
    title('raw');
    set(gca,'box','on','linewidth',1,'fontsize',10,...
        'YDir','normal','Xdir','normal');             

    subplot(142);
    errorbar(count_center(:,rr),count_sigma(:,rr),'ko',...
        'markerfacecolor','k');
    xlabel('run #');
    ylabel('estimate peak fluorescence')
    ylim([0 max(count_center(:,rr))*1.2]);

    % Spatial Compensation
    subplot(143);
    pHistB1 = bar(centers(iL),spatialN(iL),'linestyle','none',...
        'facecolor','k');
    xlim([0 max(edges)]);    
    ylabel('occurences');
    xlabel('counts per lattice site');
    hold on
    ylim([0 rawN(2)]);

    % High Counts
    yyaxis right
    pHistB2 = bar(centers(iH),spatialN(iH),'linestyle','none',...
    'FaceColor',[0.6 0 0.5]);
    set(gca,'box','on','YColor',[0.6 0 0.5]*.8);
    ylabel('occurences');
    title('with spatial compensation');
    set(gca,'box','on','linewidth',1,'fontsize',10,...
        'YDir','normal','Xdir','normal');          
   descStr=['$N_0 = ' num2str(round(mean(count_center(:,rr)))) ...
       '\pm' num2str(round(std(count_center(:,rr)))) '$'];
     text(.01,.99,descStr,'units','normalized',...
          'interpreter','latex','horizontalalignment','left',...
          'verticalalignment','cap');
      
    % Spatial Compensation with Normalization
    ax4=subplot(144);
    pHistB1 = bar(centers_norm(iL_norm),normN(iL_norm),'linestyle','none',...
        'facecolor','k');
    xlim([0 max(edges_norm)]);    
    ylabel('occurences');
    xlabel('scaled counts per lattice site');
    ylim([0 rawN(2)]);


    hold on
    % High Counts
    yyaxis right
    pHistB2 = bar(centers_norm(iH_norm),normN(iH_norm),'linestyle','none',...
    'FaceColor',[0.6 0 0.5]);
    set(gca,'box','on','YColor',[0.6 0 0.5]*.8);
    ylabel('occurences');
    title('with spatial compensation + normalized');
    set(gca,'box','on','linewidth',1,'fontsize',10,...
        'YDir','normal','Xdir','normal');
    
    i_larger_than_one = [centers_norm>=1];
    
    lvls = normN(i_larger_than_one);
    
    Ntot = 2*sum(lvls);
    mu = out(rr).PDF1_Center;
    s = out(rr).PDF1_Radius;
    foo = @(x) 1/sqrt(2*pi*s^2)*exp(-(x-mu).^2/(2*s^2));

    y=foo(centers_norm);
    
    yN= max(normN(iH_norm))*y/max(y);
    hold on
    plot(centers_norm, yN,'k--','parent',ax4,'linewidth',2);   
    
    
      if isfield(bindata,'SourceDirectory') 
        tFig=uicontrol('style','text','string',bindata(1).SourceDirectory,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left','parent',figs(rr));
        tFig.Position(4)=15;
        tFig.Position(3)=figs(rr).Position(3);
        tFig.Position(1:2)=[5 figs(rr).Position(4)-tFig.Position(4)];
      end    
      
    descStr=['$N_0 = ' num2str(round(mu,2)) ', ' ...
        '\sigma = ' num2str(round(s,2)) '$'];
      text(.01,.99,descStr,'units','normalized',...
          'interpreter','latex','horizontalalignment','left',...
          'verticalalignment','cap');
      
end


end

