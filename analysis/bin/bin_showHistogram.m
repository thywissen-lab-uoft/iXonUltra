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

if ~isfield(opts,'useManualThresh')
    opts.useManualThresh = false';
end

if ~isfield(opts,'BinSource')
    opts.BinSource = 'Zbin';
    % opts.BinSource = 'ZbinRaw';
    % opts.BinSource = 'ZbinNormalized';
end

%% Get ROI

n1 = bindata(1).LatticeBin(1).n1;
n2 = bindata(1).LatticeBin(1).n2;

if isequal(opts.ROI,'max')
   opts.ROI = [min(n1) max(n1) min(n2) max(n2)]; 
end    
R = opts.ROI;
%% Prepare Data
Zall = zeros(length(n2),length(n1),length(bindata));
for nn = 1:length(bindata)        
    Zthis = bindata(nn).LatticeBin(1).(opts.BinSource);    
    Zall(:,:,nn) =  Zthis;
    threshes(nn)=bindata.LatticeBin(1).ClusterThreshold;    
    centerval(nn) = bindata.LatticeBin(1).PDF1_Center;
end

in1i = find(n1==R(1),1);in1f = find(n1==R(2),1);    
in2i = find(n2==R(3),1);in2f = find(n2==R(4),1);
 
z=Zall;
z(z==0)=[];
[N,edges] = histcounts(z,opts.Bins);  
centers = (edges(1:end-1) + edges(2:end))/2;   

thresh = median(threshes);

centerval = median(centerval);
if isequal(opts.BinSource,'ZbinNormalized')
    thresh = thresh/centerval;
end
    iL = centers<=thresh;
iH = ~iL; 
%% Initialize Figure
if ~isfield(opts,'Parent')
    opts.Parent = figure('color','w',[100 100 1300 400],...
        'Name','BinHistogram');
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
        'w','horizontalalignment','left');
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
axImg.Position=[ax1.Position(1)+0.01 ...
    ax1.Position(2)+ax1.Position(4)-0.35*ax1.Position(4)-.01 ...
    ax1.Position(3)*0.35 ...
    ax1.Position(4)*0.35];
hImg = imagesc(n1,n2,sum(Zall,3),'parent',axImg);   
set(axImg,'visible','off','ydir','normal')
axis(axImg,'equal');
axis(axImg,'tight');
ca = [0 0 0];       
cb = [0.7 .1 .6];
cc = [linspace(ca(1),cb(1),1000)' ...
    linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
colormap(axImg,cc);

rectangle('Position',[R(1) R(3) R(2)-R(1) R(4)-R(3)],...
    'EdgeColor','r','parent',axImg)

    
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