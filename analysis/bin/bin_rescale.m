function  [bindata,hF,hF2] = bin_rescale(bindata,opts)
% This function accounts for 


opts.do_SpatialInhomogeneity =1;
opts.SpatialInhomogeneityRadius = 50; % gaussian radisu in sites
opts.SpatialInhomogeneityMaxScale = 1.5;

opts.do_Intensity = 1;
s=opts.SpatialInhomogeneityRadius;

%% Get ROI
    
    n1 = bindata(1).LatticeBin(1).n1;
    n2 = bindata(1).LatticeBin(1).n2;
    
    if isequal(opts.ROI,'max')
       opts.ROI = [min(n1) max(n1) min(n2) max(n2)]; 
    end    
    R = opts.ROI;
    in1i = find(n1==R(1),1);in1f = find(n1==R(2),1);    
    in2i = find(n2==R(3),1);in2f = find(n2==R(4),1);
 
%% Gather Data

ZbinScaled = zeros(numel(n2),numel(n1),length(bindata));

for kk=1:length(bindata)
    ZbinScaled(:,:,kk) = bindata(kk).LatticeBin(1).Zbin;
end

% Asign to ROI
ZbinScaled = ZbinScaled(in2i:in2f,in1i:in1f,:);

Zbin = ZbinScaled;

n1 = n1(in1i:in1f);
n2 = n2(in2i:in2f);
[nn1,nn2]=meshgrid(n1,n2);

ZbinScaledCopy = ZbinScaled;
ZbinScaledCopy(isnan(ZbinScaledCopy))=0;
%% Find Center of fluorescence
n1c=sum(sum(ZbinScaledCopy,3).*nn1,'all')/sum(ZbinScaledCopy,'all');
n2c=sum(sum(ZbinScaledCopy,3).*nn2,'all')/sum(ZbinScaledCopy,'all');

%% Rescale Function
gauss_map = exp(-(nn1-n1c).^2/(2*s^2)).*exp(-(nn2-n2c).^2/(2*s^2));
gauss_map_inv = 1./gauss_map;
gauss_map_inv(gauss_map_inv>opts.SpatialInhomogeneityMaxScale)=opts.SpatialInhomogeneityMaxScale;

%% Spatial Data

if opts.do_SpatialInhomogeneity
    for kk=1:size(ZbinScaled,3)
        ZbinScaled(:,:,kk) = ZbinScaled(:,:,kk).*gauss_map_inv;    
    end
    bindata(kk).ScaleRadius = s;  
end

%% Rescale Fluoresence
if opts.do_Intensity
        
    % Find x2 clusters (occupied, unoccupied)
    cluster_number=2;
    
    for kk=1:size(ZbinScaled,3)
        data=ZbinScaled(:,:,kk);
        data=data(:);
        data(isnan(data))=[];
        data(data==0)=[];
        [idx,c,sumD,D]=kmeans(data,cluster_number);

    % Sort by centroid
    [c,inds]=sort(c,'ascend');
    sumD=sumD(inds);

    dd = zeros(numel(c),1);
    thresh = zeros(numel(c),1);
        nCluster=zeros(numel(c),1);

        for bb=1:numel(c)
            ind = inds(bb);
            nThis = sum(idx==ind);
            dd(bb) = sqrt(sumD(bb)/nThis);

            thresh(bb) = round(max(data(idx==ind)));
            nCluster(bb)=nThis;
        end
          
        
        z2 = data(data>thresh(1));
        % The n=1 (atom) distribution is fit with a simple gaussian
        % distribution
        % pd1 = fitdist(z2,'normal'); % This works too, but use MLE to make it
        % the same output
        [pdf1_c,pdf1_cint] = mle(z2,'distribution','normal');
        
        I0=pdf1_c(1); % center and then the 
        
                
        ZbinScaled(:,:,kk)=ZbinScaled(:,:,kk)/I0;
        
        out.ClusterNumber = cluster_number;
        out.ClusterThreshold = thresh;
        out.ClusterCentroid = c;
        out.ClusterRadius = dd;        
        out.ClusterScale = I0;
        
        bindata(kk).ScaledCentroid = I0;
        bindata(kk).ScaledCentroidRadius = pdf1_c(2)/I0;
        bindata(kk).ScaledThreshold = (I0-2*pdf1_c(2))/I0;
    end      
end

Zscaledplot = ZbinScaled;
ZscaledplotBar = mean(Zscaledplot,3);
Zscaledplot(Zscaledplot<1e-2)=[];
Zscaledplot(isnan(Zscaledplot))=[];


scaledEdges = linspace(0,3,50);

scaledN = histcounts(Zscaledplot,scaledEdges); 
scaledCenters = (scaledEdges(1:end-1) + scaledEdges(2:end))/2;      

scaledThreshold = mean([bindata.ScaledThreshold]);

scaled_iL = scaledCenters<=scaledThreshold;
scaled_iH = ~scaled_iL;

I0_bar = mean([bindata.ScaledCentroid]);
I0_std = std([bindata.ScaledCentroid]);

I0_radius = mean([bindata.ScaledCentroidRadius]);
I0_radius_std = std([bindata.ScaledCentroidRadius]);

descStr=['$N_0 = ' num2str(round(I0_bar)) '\pm' num2str(round(I0_std)) '$' ...
    newline '$\sigma_{\mathrm{counts}} = ' num2str(round(100*I0_radius)) '\% \pm' num2str(round(100*I0_radius_std)) '\% $' ...
    newline '$\mathrm{gauss~radius:}' num2str(s) '~\mathrm{sites}$'];

%% Analysis on Raw Data

Zbinplot = Zbin;
ZbinplotBar = mean(Zbinplot,3);
Zbinplot(Zbinplot<20)=[];
Zbinplot(isnan(Zbinplot))=[];
[idx,c,sumD,D]=kmeans(Zbinplot(:),cluster_number);
% Sort by centroid
[c,inds]=sort(c,'ascend');
sumD=sumD(inds);

dd = zeros(numel(c),1);
thresh = zeros(numel(c),1);
nCluster=zeros(numel(c),1);

for bb=1:numel(c)
    ind = inds(bb);
    nThis = sum(idx==ind);
    dd(bb) = sqrt(sumD(bb)/nThis);

    thresh(bb) = round(max(Zbinplot(idx==ind)));
    nCluster(bb)=nThis;
end

rawEdges = linspace(0,3*c(2),50);

rawN = histcounts(Zbinplot,rawEdges); 
centers = (rawEdges(1:end-1) + rawEdges(2:end))/2;      

iL = centers<=thresh(1);
iH = ~iL; 

%%
hF = figure;
hF.Color='w';
hF.Position=[50 50 900 600];

 

ca = [0 0 0];       
    cb = [0.7 .1 .6];
    cc = [linspace(ca(1),cb(1),1000)' ...
        linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
    colormap(hF,cc);
    hold on


subplot(221);
imagesc(n1,n2,ZbinplotBar);
xlabel('n1 sites');
ylabel('n2 sites');
title(['uncompensated avg image; N=' num2str(length(bindata))]);
set(gca,'box','on','linewidth',1,'fontsize',10,...
    'YDir','normal','Xdir','normal');
colorbar;
caxis([0 c(2)]);

subplot(222);
pHistB1 = bar(centers(iL),rawN(iL),'linestyle','none',...
    'facecolor','k');
xlim([0 max(rawEdges)]);    
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
imagesc(n1,n2,ZscaledplotBar);
xlabel('n1 sites');
ylabel('n2 sites');
title(['compensated avg image; N=' num2str(length(bindata))]);
set(gca,'box','on','linewidth',1,'fontsize',10,...
    'YDir','normal','Xdir','normal');
colorbar;
caxis([0 1]);

subplot(224);
pHistB1 = bar(scaledCenters(scaled_iL),scaledN(scaled_iL),'linestyle','none',...
    'facecolor','k');
xlim([0 max(scaledEdges)]);    
ylabel('occurences');
xlabel('scaled counts per lattice site');
hold on
% High Counts
yyaxis right
pHistB2 = bar(scaledCenters(scaled_iH),scaledN(scaled_iH),'linestyle','none',...
'FaceColor',[0.6 0 0.5]);
set(gca,'box','on','YColor',[0.6 0 0.5]*.8);
ylabel('occurences');
title('compensated histogram');
set(gca,'box','on','linewidth',1,'fontsize',10,...
    'YDir','normal','Xdir','normal');

      if isfield(bindata,'SourceDirectory') 
        tFig=uicontrol('style','text','string',bindata(1).SourceDirectory,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left','parent',hF);
        tFig.Position(4)=15;
        tFig.Position(3)=hF.Position(3);
        tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
      end    
    
    
      text(.95,.95,descStr,'units','normalized',...
          'interpreter','latex','horizontalalignment','right',...
          'verticalalignment','cap');
      
%% Other Figure

      P=[bindata.Params];
    X = [P.(opts.xVar)];
    
    Y =  [bindata.ScaledCentroid];
    Ye = [bindata.ScaledCentroidRadius].*[bindata.ScaledCentroid];
   
    
hF2 = figure;
hF2.Color='w';
hF2.Position=[1000 50 600 400];


    ax =axes;   
    
    co=get(gca,'colororder');

    hold on

  
    errorbar(X,Y,2*Ye,'o','linewidth',2,'markersize',10,'markerfacecolor',co(1,:),...
        'markeredgecolor',co(1,:)*.5);
    xlabel(opts.xVar)
    if isequal(opts.xVar,'ExecutionDate')
        datetick x
    end
    ylabel('fluorescence/atom');
    set(gca,'box','on','linewidth',1,'fontsize',10);
    grid on;
    yL=get(gca,'YLim');
    set(gca,'YLim',[0 yL(2)]);
    legend({'$\mathrm{compensated~centroid} \pm 2\sigma$'},'interpreter','latex');
%     
%              bindata(kk).ScaledCentroidRadius = pdf1_c(2)/I0;
%         bindata(kk).ScaledThreshold = (I0-2*pdf1_c(2))/I0;
%% Assign Outputs

for kk=1:length(bindata)
   bindata(kk).LatticeBin(1).ZbinScaled = ZbinScaled(:,:,kk);
   bindata(kk).LatticeBin(1).ScaledThreshold = scaledThreshold;
end

end

