function hF=ixon_showRawCountHistogram(ixondata,xVar,opts)
global ixon_imgdir
imgnum=opts.ImageNumber;

params=[ixondata.Params];
acq=[ixondata.AcquisitionInformation];

xData=1:length(ixondata);
if ismember(xVar,fieldnames(params))
    xData=[params.(xVar)];
end

if ismember(xVar,fieldnames(acq))
    xData=[acq.(xVar)];
end



%% Finding Global limits

gMin=inf;
gMax=-inf;
for kk=1:length(ixondata)
    img=ixondata(kk).RawImages(:,:,imgnum);
    vals=sort(img(:),'descend');
    
    gMin=min([gMin vals(end-opts.Outliers(1))]);
    gMax=max([gMax vals(opts.Outliers(2))]);

%     gMin=min([gMin ixondata(kk).Raw(imgnum).Minimum]);
%     gMax=max([gMax ixondata(kk).Raw(imgnum).Maximum]);
end

    
%% Make Fgiure

strs=strsplit(ixon_imgdir,filesep);
str=[strs{end-1} filesep strs{end}];



hF=figure('Name', ['Raw Histogram : ' str], 'Visible', 'On', ...
    'NumberTitle','off','color','w','MenuBar','none','units','pixels',...
    'Resize','off');  
hF.Position(1)=20;
hF.Position(2)=50;
hF.Position(3)=1850;
hF.Position(4)=1000;
clf;

t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];


for ii=1:length(ixondata)
    % Create axes object
    thisAxes=axes('parent',hF,'units','pixels');
    set(thisAxes,'FontSize',8,'XMinorTick','on','YMinorTick','on',...
        'Box','on','YGrid','on','XGrid','on','units','pixels',...
        'YAxisLocation','right');
    hold on;     
    cla;
    [a,b,c,d]=getAxesPos(ii,length(ixondata),...
        hF.Position(3),hF.Position(4));
    thisAxes.Position(1)=a;
    thisAxes.Position(2)=b;
    thisAxes.Position(3)=c;
    thisAxes.Position(4)=d;   


    img=ixondata(ii).RawImages(:,:,imgnum);
    vals=sort(img(:),'descend');
    if opts.GlobalLimits
        histogram(img,'BinWidth',opts.BinWidth,'BinLimits',[gMin gMax],...
            'edgealpha',0.0); 
        xlim([gMin-10 gMax+10]);
    else
        histogram(img,'BinWidth',opts.BinWidth); 
        set(gca,'XLim',[vals(end-opts.Outliers(1)) vals(opts.Outliers(2))],...
            'edgealpha',0.25);
    end



    str=['{\bf sum: }' ...
        num2str(ixondata(ii).Raw(imgnum).Total,'%.2e') ...
        '{\bf median: }' ...
        num2str(round(ixondata(ii).Raw(imgnum).Median,1)) ...
        '{\bf \sigma: }'  ...
        num2str(round(ixondata(ii).Raw(imgnum).Std,1))];    
    
        
    % Plot the data
    try  
        % Draw the analysis string box
        text(thisAxes.Position(3)-1, thisAxes.Position(4), str, 'Units', 'pixels',...
            'FontSize', 8,...
            'verticalalignment','cap','horizontalalignment','right',...
            'backgroundcolor',[1 1 1 .5]); 

        % Draw the iteration number and variable value
         tt=text(3, thisAxes.Position(4)-1, ...
             ['{\bf(' num2str(ii) ')' newline ...
             num2str(xData(ii)) '}'], ...
             'Units', 'pixels',...
             'FontSize', 8,...
             'verticalalignment','cap','HorizontalAlignment','left',...
             'backgroundcolor',[1 1 1 .5]);     
    end

    set(gca,'YScale',opts.YScale);
    drawnow

    
end      
disp('Done!');
end

function [axX,axY,axWidth,axHeight]=getAxesPos(nInd,nTot,xSize,ySize)
nInd=nInd-1;
yTop=30;
yBot=30;

xLeft=10;
xRight=20;

ySpace=25;
xSpace=20;

nRow=ceil(sqrt(nTot));

axHeight=(ySize-yTop-yBot-ySpace*(nRow-1))/nRow;
axWidth=(xSize-xLeft-xRight-xSpace*(nRow-1))/nRow;

axX=xLeft+(axWidth+xSpace)*mod(nInd,nRow);
axY=(ySize-yTop-axHeight)-floor(nInd/nRow)*(axHeight+ySpace);
end
