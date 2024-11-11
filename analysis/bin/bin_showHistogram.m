function  fig = bin_showHistogram(bindata,opts)
    
%% Options
if nargin == 1
    opts = struct;
end

if ~isfield(opts,'ROI')
    opts.ROI = 'max';
end  

if ~isfield(opts,'Bins')
    opts.Bins = 80;
end   

if ~isfield(opts,'doAnimate')
   opts.doAnimate = 0; 
end

if ~isfield(opts,'ManualThresh')
   opts.ManualThresh = 2000; 
end

if ~isfield(opts,'ImageNum')
    opts.ImageNum = 1;
end

if ~isfield(opts,'useManualThresh')
    opts.useManualThresh = false';
end

if ~isfield(opts,'BinSource')
    opts.BinSource = 'Zbin';
end

%% Colormap
ca = [1 1 1];
cb = [0.6 0 .5];
white2purple = [linspace(ca(1),cb(1),1000)' ...
    linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];


ca = [0 0 0];       
cb = [0.7 .1 .6];
black2purple = [linspace(ca(1),cb(1),1000)' ...
    linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
%% Get ROI

n1 = bindata(1).LatticeBin(opts.ImageNum).n1;
n2 = bindata(1).LatticeBin(opts.ImageNum).n2;

[nn1,nn2]=meshgrid(n1,n2);


%% Acrew all data
Zall = zeros(length(n2),length(n1),length(bindata));
for nn = 1:length(bindata)        
    Zthis = bindata(nn).LatticeBin(opts.ImageNum).(opts.BinSource);    
    Zall(:,:,nn) =  Zthis;
    threshes(nn)=bindata.LatticeBin(opts.ImageNum).ClusterThreshold;    
    centerval(nn) = bindata.LatticeBin(1).PDF1_Center;

    Zthis(isnan(Zthis))=0;

    n1c(nn) = sum(Zthis.*nn1,'all')/sum(Zthis,'all');
    n2c(nn) = sum(Zthis.*nn2,'all')/sum(Zthis,'all');

    n1c_sq(nn) = sum(Zthis.*nn1.^2,'all')/sum(Zthis,'all');
    n2c_sq(nn) = sum(Zthis.*nn2.^2,'all')/sum(Zthis,'all');

    n1_sigma(nn) = sqrt(n1c_sq-n1c^2);
    n2_sigma(nn) = sqrt(n2c_sq-n2c^2);
end

n_sigma_med =[median(n1_sigma) median(n2_sigma)];
nc_med = [median(n1c) median(n2c)];

n1_lim = nc_med(1)+2*[-1 1]*n_sigma_med(1);
n2_lim = nc_med(2)+2*[-1 1]*n_sigma_med(2);


z=Zall;
z(z==0)=[];
[N,edges] = histcounts(z,opts.Bins);  
centers = (edges(1:end-1) + edges(2:end))/2;   

thresh = median(threshes);
centerval = median(centerval);
iL = centers<=thresh;
iH = ~iL; 
Nall={};

for nn=1:length(bindata)
    zthis = Zall(:,:,nn); zthis=zthis(:);
    zthis(zthis==0)=[];
    zthis(isnan(zthis))=[];
    Nall{nn}=histcounts(z,edges);
end

%% Initialize Figure
if ~isfield(opts,'Parent')
    opts.Parent = figure('color','w','Position',[100 100 550 400],...
        'Name','BinHistogram','NumberTitle','off');
    fig = opts.Parent;
else
    fig = opts.Parent;
    for kk=1:length(fig.Children)
        delete(fig.Children(1))
    end
end
    
if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
    tFig=uicontrol('style','text','string',opts.FigLabel,...
        'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left','parent',fig);
    tFig.Position(4)=tFig.Extent(4);
    tFig.Position(3)=400;
    tFig.Position(1:2)=[1 1];
end    

% Histogram Axis
ax1 = axes('parent',opts.Parent);

% Low Counts
pHistB1 = bar(centers(iL),N(iL),'linestyle','none',...
    'facecolor','k','parent',ax1);
xlim(ax1,[0 max(edges)]);    
ylabel(ax1,'occurences');
xlabel(ax1,'counts per lattice site');
    
% High Counts
yyaxis(ax1,'right');
pHistB2 = bar(centers(iH),N(iH),'linestyle','none',...
    'FaceColor',[0.6 0 0.5],'parent',ax1);
ylabel(ax1,'occurences');
set(ax1,'box','on','YColor',[0.6 0 0.5]*.8);
xlim(ax1,[0 max(edges)]);    

% Image Axis
axImg = axes('parent',opts.Parent);
axImg.Position=[ax1.Position(1) ...
    ax1.Position(2)+ax1.Position(4)-0.4*ax1.Position(4) ...
    ax1.Position(3)*0.4 ...
    ax1.Position(4)*0.4];
hImg = imagesc(n1,n2,sum(Zall,3),'parent',axImg);   
set(axImg,'visible','off','ydir','normal')
axis(axImg,'equal');
axis(axImg,'tight');
colormap(axImg,black2purple);
% caxis(axImg,[0 thresh]);
set(axImg,'XLim',n1_lim,'YLim',n2_lim);
colorbar(axImg);
% rectangle('Position',[R(1) R(3) R(2)-R(1) R(4)-R(3)],...
%     'EdgeColor','r','parent',axImg)

    
    %% Populate with data
    
    % if opts.doAnimate
    %     for nn = 1:length(bindata)
    %         [N,~] = histcounts(Zall(in2i:in2f,in1i:in1f,nn),opts.Bins);  
    %         set(pHistB1,'Xdata',centers(iL),'Ydata',N(iL));
    %         set(pHistB2,'Xdata',centers(iH),'Ydata',N(iH));
    %         set(hImg,'Cdata',Zall(:,:,nn));
    %            str1 = ['ROI : [' num2str(R(1)) ' ' num2str(R(2)) ' ' ...
    %     num2str(R(3)) ' ' num2str(R(4)) ']' ...
    %     newline 'image ' num2str(nn) '/' num2str(length(bindata))];
    % 
    %         set(t1,'String',str1);
    % 
    %         frame=getframe(hF);
    %         im = frame2im(frame);
    %         [A,map] = rgb2ind(im,256);              
    %         filename = fullfile(opts.saveDir,opts.filename);            
    %         switch nn
    %             case 1
    %                 imwrite(A,map,filename,'gif','LoopCount',...
    %                     Inf,'DelayTime',1);
    %             case length(bindata)
    %                 imwrite(A,map,filename,'gif','WriteMode',...
    %                     'append','DelayTime',1);
    %             otherwise
    %                 imwrite(A,map,filename,'gif','WriteMode',...
    %                     'append','DelayTime',.1);
    %         end
    %     end
    % end
    

end