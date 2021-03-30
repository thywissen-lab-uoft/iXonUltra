function hF=showRawCountHistogram(atomdata,xVar,opts)
global imgdir
imgnum=opts.ImageNumber;

params=[atomdata.Params];
acq=[atomdata.AcquisitionInformation];

% keyboard
xData=1:length(atomdata);
if ismember(xVar,fieldnames(params))
    xData=[params.(xVar)];
end

if ismember(xVar,fieldnames(acq))
    xData=[acq.(xVar)];
end



%% Finding Global limits

gMin=inf;
gMax=-inf;
for kk=1:length(atomdata)
    gMin=min([gMin atomdata(kk).Raw(imgnum).Minimum]);
    gMax=max([gMax atomdata(kk).Raw(imgnum).Maximum]);
end

    
%% Make Fgiure

strs=strsplit(imgdir,filesep);
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


for ii=1:length(atomdata)
    % Create axes object
    thisAxes=axes('parent',hF,'units','pixels');
    set(thisAxes,'FontSize',8,'XMinorTick','on','YMinorTick','on',...
        'Box','on','YGrid','on','XGrid','on','units','pixels',...
        'YTickLabel',{});
    hold on;     
    cla;
    [a,b,c,d]=getAxesPos(ii,length(atomdata),...
        hF.Position(3),hF.Position(4));
    thisAxes.Position(1)=a;
    thisAxes.Position(2)=b;
    thisAxes.Position(3)=c;
    thisAxes.Position(4)=d;   
    
    
    if opts.GlobalLimits
        histogram(atomdata(ii).RawImages(:,:,imgnum),'BinWidth',opts.BinWidth,'BinLimits',[gMin gMax]); 
    else
        histogram(atomdata(ii).RawImages(:,:,imgnum),'BinWidth',opts.BinWidth); 
    end

   

    str=['{\bf sum: }' ...
        num2str(atomdata(ii).Raw(imgnum).Total,'%.2e') ...
        '{\bf median: }' ...
        num2str(round(atomdata(ii).Raw(imgnum).Median,1)) ...
        '{\bf \sigma: }'  ...
        num2str(round(atomdata(ii).Raw(imgnum).Std,1))];    
    
        
    % Plot the data
    try

    
    

    % Draw the analysis string box
    text(thisAxes.Position(3)-1, thisAxes.Position(4), str, 'Units', 'pixels',...
        'FontSize', 8,...
        'verticalalignment','cap','horizontalalignment','right'); 

    % Draw the iteration number and variable value
     text(3, thisAxes.Position(4)-1, ...
         ['{\bf(' num2str(ii) ')' newline ...
         num2str(xData(ii)) '}'], ...
         'Units', 'pixels',...
         'FontSize', 8,...
         'verticalalignment','cap','HorizontalAlignment','left'); 
    
    end
    set(gca,'YScale',opts.YScale);

end      
disp('Done!');
end

function [axX,axY,axWidth,axHeight]=getAxesPos(nInd,nTot,xSize,ySize)
nInd=nInd-1;
yTop=30;
yBot=30;

xLeft=20;
xRight=20;

ySpace=25;
xSpace=10;

nRow=ceil(sqrt(nTot));

axHeight=(ySize-yTop-yBot-ySpace*(nRow-1))/nRow;
axWidth=(xSize-xLeft-xRight-xSpace*(nRow-1))/nRow;

axX=xLeft+(axWidth+xSpace)*mod(nInd,nRow);
axY=(ySize-yTop-axHeight)-floor(nInd/nRow)*(axHeight+ySpace);
end
