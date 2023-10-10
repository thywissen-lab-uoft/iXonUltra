function hFs=ixon_showProfile(ixondata,direction,xVar,opts)
    
pMax=36;

rNum = opts.ROINum;
style = opts.Style;

if nargin == 4 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end
    
%% Make Fgiure

clear hFs
for kk=1:(ceil(length(ixondata)/pMax))    
    nStart=(kk-1)*pMax+1;
    nEnd=min([pMax*kk length(ixondata)]);
    fprintf(['Showing OD profile ' direction ' ROI ' ...
        num2str(rNum) ' ' num2str(nStart) ' to ' num2str(nEnd) ' ... ']);

    ixondataSUB=ixondata(nStart:nEnd);  
    
    hFs(kk)=figure('Name', [pad(['Cut ' direction ' R' num2str(rNum) ' ' num2str(kk)],20) FigLabel], 'Visible', 'On', ...
        'NumberTitle','off','color','w','MenuBar','none','units','pixels',...
        'Resize','off'); 
    hF=hFs(kk);
    hF.Position(1)=20;
    hF.Position(2)=50;
    hF.Position(3)=1850;
    hF.Position(4)=1000;
    clf;

    t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    t.Position(4)=t.Extent(4);
    t.Position(3)=hF.Position(3);
    t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];


    for ii=1:length(ixondataSUB)
        % Create axes object
        ax=axes('parent',hF,'units','pixels');
        set(ax,'FontSize',8,'XMinorTick','on','YMinorTick','on',...
            'Box','on','YGrid','on','XGrid','on','units','pixels',...
            'YTickLabel',{});
        hold on;     
        cla;
        [a,b,c,d]=getAxesPos(ii,length(ixondataSUB),...
            hF.Position(3),hF.Position(4));
        ax.Position = [a b c d];
        
        % Get this ROI
        ROI=ixondataSUB(ii).ROI(rNum,:);

        
        ix_1 = find(ixondataSUB(ii).X>=ROI(1),1);
        ix_2 = find(ixondataSUB(ii).X>=ROI(2),1);
        iy_1 = find(ixondataSUB(ii).Y>=ROI(3),1);
        iy_2 = find(ixondataSUB(ii).Y>=ROI(4),1);            
        x = ixondataSUB(ii).X(ix_1:ix_2);
        y = ixondataSUB(ii).Y(iy_1:iy_2);   
        z = ixondataSUB(ii).Z(iy_1:iy_2,ix_1:ix_2,rNum);         
        
        
        % Get the data
%         x=ixondataSUB(ii).X;y=ixondataSUB(ii).Y;  
%         z=ixondataSUB(ii).OD; % N counts

        % Get data over the selected ROI
%         x=x(ROI(1):ROI(2));y=y(ROI(3):ROI(4));
%         z=z(ROI(3):ROI(4),ROI(1):ROI(2));
        
        switch direction
            case 'X'
                X=x;
            case 'Y'
                X=y;
            otherwise
                error('invalid plot direction');                               
        end
        
        % Mesh grid for fits
        [xx,yy]=meshgrid(x,y);
        
        Yc = [];
        Xc = [];
        
        % Get the gaussian fit
        clear gaussFit
        doGauss = 0;
        if isfield(ixondataSUB(ii),'GaussFit') && ~isfield(ixondataSUB(ii),'FermiFit')
            gaussFit  = ixondataSUB(ii).GaussFit{rNum};
            Yc(end+1) = gaussFit.Yc;
            Xc(end+1) = gaussFit.Xc;
            doGauss = 1;
        end    
        
        % Get the box count
        clear doBox
        doBox = 0;        
        if isfield(ixondataSUB(ii),'BoxCount') && ~(doGauss)
            Yc(end+1) = ixondataSUB(ii).BoxCount(rNum).Yc;
            Xc(end+1) = ixondataSUB(ii).BoxCount(rNum).Xc;
            
            if Yc<y(1) || Yc>y(end)
                Yc = mean(y);
            end
            
            if Xc<x(1) || Xc>x(end)
                Xc = mean(x);
            end
            
            doBox = 1;
        end   
        
        % Find index to plot against
        Yc = mean(Yc);        
        iY = find(y>=Yc,1);   
        
        Xc = mean(Xc);
        iX = find(x>=Xc,1);
        
        %%%% Get gauss profile %%%%
        if doGauss   
            zzF_gauss = feval(gaussFit,xx,yy);

            if isequal(direction,'X') && isequal(style,'cut')
                YF_gauss = zzF_gauss(iY,:);
            end
            
            if isequal(direction,'X') && isequal(style,'sum')
                YF_gauss = sum(zzF_gauss,1);
            end
            
            if isequal(direction,'Y') && isequal(style,'cut')
                YF_gauss = zzF_gauss(:,iX);
            end
            
            if isequal(direction,'Y') && isequal(style,'sum')
                YF_gauss = sum(zzF_gauss,2);
            end            
        end
        
        
        %%%% Get the data %%%%
        if isequal(direction,'X') && isequal(style,'cut')
            Y_data = z(iY,:);
        end

        if isequal(direction,'X') && isequal(style,'sum')
            Y_data = sum(z,1);
        end

        if isequal(direction,'Y') && isequal(style,'cut')
            Y_data = z(:,iX);
        end

        if isequal(direction,'Y') && isequal(style,'sum')
            Y_data = sum(z,2);
        end  
        %%%% Plot the fits %%%%
        if doGauss
            plot(X,YF_gauss,'r','LineWidth',2);
        end
        
        
        % Plot the data
        plot(X,Y_data,'k-');

        % Adjust limits
        xlim([X(1) X(end)]);   
        ax.YLim(1)=min([0 min(Y_data)]);
        ax.YLim(2)=max([max(Y_data)*1.5 0]);   
        
        % Draw the analysis string box
        iterNum=(kk-1)*pMax+ii;
        
        % x variable string
        if isequal(xVar,'ExecutionDate')
            xstr = datestr(ixondataSUB(ii).Params.(xVar),'mm/DD HH:MM:ss');
        else
            xstr = num2str(ixondataSUB(ii).Params.(xVar));
        end

        % Draw the iteration number and variable value
        text(3, ax.Position(4)-2, ...
            ['{\bf(' num2str(iterNum) ')' newline ...
            xstr '}'], ...
            'Units', 'pixels',...
            'FontSize', 8,...
            'verticalalignment','cap','HorizontalAlignment','left');
        
        %%%% Determine String Label %%%%
        lstr = [direction ' ' style];
        if isequal(style,'cut') && isequal(direction,'X')
           lstr = [lstr ' @ y = ' num2str(round(Yc)) ' px']; 
        end
        
        if isequal(style,'cut') && isequal(direction,'Y')
           lstr = [lstr ' @ x = ' num2str(round(Xc)) ' px']; 
        end
        
        if doGauss
            if isequal(direction,'X')
               lstr = [lstr newline ...
                   'gauss {\bf c: }'  num2str(round(gaussFit.Xc)) ...
                    'px  ' ...
                    '{\bf \sigma: }' num2str(round(gaussFit.Xs)) ...
                    'px'];   
            else
                lstr = [lstr newline ...
                   'gauss {\bf c: }'  num2str(round(gaussFit.Yc)) ...
                    'px  ' ...
                    '{\bf \sigma: }' num2str(round(gaussFit.Ys)) ...
                    'px'];   
            end
        end
        

        % Draw the analysis string box
        text(ax.Position(3)-1, ax.Position(4)-2, lstr, 'Units', 'pixels',...
            'FontSize', 8,...
            'verticalalignment','cap','horizontalalignment','right'); 
        
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
