function hFs=ixon_showGaussProfile(ixondata,direction,style,xVar)
global ixon_imgdir

pMax=36;


switch nargin
    case 1
        direction='X';
        xVar='timestamp';
        style='cut';
        
    case 2
        xVar='timestamp';
        style='cut';
    case 3
        xVar='timestamp';

end
    
%% Make Fgiure


clear hFs
for kk=1:(ceil(length(ixondata)/pMax))
    
    nStart=(kk-1)*pMax+1;
    nEnd=min([pMax*kk length(ixondata)]);
    fprintf(['Showing gauss profile ' direction ' ' num2str(nStart) ' to ' num2str(nEnd) ' ... ']);

    ixondataSUB=ixondata(nStart:nEnd);
    
    
    strs=strsplit(ixon_imgdir,filesep);
    str=[strs{end-1} filesep strs{end}];
    
    hFs(kk)=figure('Name', [pad(['Ixon Cut ' direction ' ' num2str(kk)],20) str], 'Visible', 'On', ...
        'NumberTitle','off','color','w','MenuBar','none','units','pixels',...
        'Resize','off'); 
    hF=hFs(kk);
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


    for ii=1:length(ixondataSUB)
        % Create axes object
        thisAxes=axes('parent',hF,'units','pixels');
        set(thisAxes,'FontSize',8,'XMinorTick','on','YMinorTick','on',...
            'Box','on','YGrid','on','XGrid','on','units','pixels',...
            'YTickLabel',{});
        hold on;     
        cla;
        [a,b,c,d]=getAxesPos(ii,length(ixondataSUB),...
            hF.Position(3),hF.Position(4));
        thisAxes.Position(1)=a;
        thisAxes.Position(2)=b;
        thisAxes.Position(3)=c;
        thisAxes.Position(4)=d;   

        ROI=ixondataSUB(ii).ROI(1,:);

        % Get the data
        x=ixondataSUB(ii).X;    % X vector
        y=ixondataSUB(ii).Y;    % Y vector
        z=ixondataSUB(ii).Z; % N counts

        % Get data over the selected ROI
        x=x(ROI(1):ROI(2));
        y=y(ROI(3):ROI(4));
        z=z(ROI(3):ROI(4),ROI(1):ROI(2));

        % Get the gaussian fit
        fout=ixondataSUB(ii).GaussFit{1};


        % Evaluvate the fit for doing numerical projection
        [xx,yy]=meshgrid(x,y);
        zzF=feval(fout,xx,yy);      

        indy=find(round(fout.Yc)==y);           % Y center
        indx=find(round(fout.Xc)==x);           % X center      

        if isequal(direction,'X')
            switch style
                case 'cut'
                    Y=z(indy,:); % Z(Xc,x)
                    YF=zzF(indy,:);
                case 'sum'
                    Y=sum(z,1); % Z(Xc,x)
                    YF=sum(zzF,1);
            end


            X=x;        



            str=['{\bf x_c: }'  num2str(round(fout.Xc)) ...
                'px'];   
        else
            switch style
                case 'cut'
                    Y=z(:,indx); % Z(Xc,x)
                    YF=zzF(:,indx);
                case 'sum'
                    Y=sum(z,2); % Z(Xc,x)
                    YF=sum(zzF,2);
            end
            X=y;        

  

            str=['{\bf y_c: }'  num2str(round(fout.Yc)) ...
                'px'];    
        end

        [cface,cedge] = ixoncolororder(1);

        % Plot the data
        try
            plot(X,YF,'-','LineWidth',3,'color',[cedge .5]);
            plot(X,Y,'k-');  

        % Set the limits
        xlim([X(1) X(end)]);   
        thisAxes.YLim(1)=min([0 min(Y)]);
        thisAxes.YLim(2)=max([max(Y)*1.5 0]);       


        % Draw the analysis string box
        text(thisAxes.Position(3)-1, thisAxes.Position(4), str, 'Units', 'pixels',...
            'FontSize', 8,...
            'verticalalignment','cap','horizontalalignment','right'); 

        iterNum=(kk-1)*pMax+ii;
        
        % Draw the iteration number and variable value
        text(3, thisAxes.Position(4)-1, ...
            ['{\bf(' num2str(iterNum) ')' newline ...
            num2str(ixondataSUB(ii).Params.(xVar)) '}'], ...
            'Units', 'pixels',...
            'FontSize', 8,...
            'verticalalignment','cap','HorizontalAlignment','left'); 
        end
    end      
    disp('done.');
    
end
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
