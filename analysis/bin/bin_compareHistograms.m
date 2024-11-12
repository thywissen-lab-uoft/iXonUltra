function  figs = bin_compareHistograms(bindata)
% This function accounts for 


%% Init Data
Npics = length(bindata(1).LatticeBin);
count_center = zeros(length(bindata),Npics);
count_sigma  = zeros(length(bindata),Npics);

%% Get all Count Centers
for kk=1:length(bindata)
    for rr=1:length(bindata(kk).LatticeBin)
        count_center(kk,rr) = bindata(kk).LatticeBin(rr).PDF1_Center;
        count_sigma(kk,rr) = bindata(kk).LatticeBin(rr).PDF1_Radius; 
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

NthreshPlot = median(count_center(:))-3*median(count_sigma(:));
edges=linspace(0,3*NthreshPlot,50);

med_norm_sigma = mean([count_sigma(:)./count_center(:)]);
std_norm_sigma = std([count_sigma(:)./count_center(:)]);

NthreshNormPlot = (1-3*med_norm_sigma);
edges_norm = linspace(0,3*NthreshNormPlot,50);


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
    figs(rr).Position=[50 50 900 600];  
    colormap(figs(rr),cc);        
    hold on
    
    zr=raw{rr}(:);
    zr(zr==0)=[];
    rawN = histcounts(zr,edges); 
    centers = (edges(1:end-1) + edges(2:end))/2;    
    iL = centers<=NthreshPlot;
    iH = ~iL; 
    
    zn=normalized{rr}(:);
    zn(zn==0)=[];
    normN = histcounts(zn,edges_norm); 
    centers_norm = (edges_norm(1:end-1) + edges_norm(2:end))/2;    
    iL_norm = centers_norm<=NthreshNormPlot;
    iH_norm = ~iL_norm; 

    subplot(221);
    imagesc(n1,n2,mean(raw{rr},3));
    xlabel('n1 sites');
    ylabel('n2 sites');
    title(['raw avg image; N=' num2str(length(bindata))]);
    set(gca,'box','on','linewidth',1,'fontsize',10,...
        'YDir','normal','Xdir','normal');
    colorbar;
    caxis([0 NthreshPlot]);

    subplot(222);
    pHistB1 = bar(centers(iL),rawN(iL),'linestyle','none',...
        'facecolor','k');
    xlim([0 max(edges)]);    
    ylabel('occurences');
    xlabel('counts per lattice site');
    hold on
    % High Counts
    yyaxis right
    pHistB2 = bar(centers(iH),rawN(iH),'linestyle','none',...
    'FaceColor',[0.6 0 0.5]);
    set(gca,'box','on','YColor',[0.6 0 0.5]*.8);
    ylabel('occurences');
    title('uncompensated histogram');
    set(gca,'box','on','linewidth',1,'fontsize',10,...
        'YDir','normal','Xdir','normal');                     

    subplot(223);
    imagesc(n1,n2,mean(normalized{rr},3));
    xlabel('n1 sites');
    ylabel('n2 sites');
    title(['compensated avg image; N=' num2str(length(bindata))]);
    set(gca,'box','on','linewidth',1,'fontsize',10,...
        'YDir','normal','Xdir','normal');
    colorbar;
    caxis([0 NthreshNormPlot]);

    subplot(224);
    pHistB1 = bar(centers_norm(iL_norm),normN(iL_norm),'linestyle','none',...
        'facecolor','k');
    xlim([0 max(edges_norm)]);    
    ylabel('occurences');
    xlabel('scaled counts per lattice site');
    hold on
    % High Counts
    yyaxis right
    pHistB2 = bar(centers_norm(iH_norm),normN(iH_norm),'linestyle','none',...
    'FaceColor',[0.6 0 0.5]);
    set(gca,'box','on','YColor',[0.6 0 0.5]*.8);
    ylabel('occurences');
    title('compensated histogram');
    set(gca,'box','on','linewidth',1,'fontsize',10,...
        'YDir','normal','Xdir','normal');

      if isfield(bindata,'SourceDirectory') 
        tFig=uicontrol('style','text','string',bindata(1).SourceDirectory,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left','parent',figs(rr));
        tFig.Position(4)=15;
        tFig.Position(3)=figs(rr).Position(3);
        tFig.Position(1:2)=[5 figs(rr).Position(4)-tFig.Position(4)];
      end    
   descStr=['$N_0 = ' num2str(round(mean(count_center(:,rr)))) '\pm' num2str(round(std(count_center(:,rr)))) '$' ...
     newline '$\sigma_{\mathrm{counts}} = ' num2str(round(100*med_norm_sigma)) '\% \pm' num2str(round(100*std_norm_sigma)) '\% $'];
 
      text(.05,.95,descStr,'units','normalized',...
          'interpreter','latex','horizontalalignment','left',...
          'verticalalignment','cap');
end


end

