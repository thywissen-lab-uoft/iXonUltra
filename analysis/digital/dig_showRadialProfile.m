function hFs=dig_showRadialProfile(digdata,opts)
    
pMax=25;

if nargin==1
    opts=struct;
end

  if ~isfield(opts,'rMaxShow')
        opts.rMaxShow = 80;
  end
  
  if ~isfield(opts,'showDevParametrization')
     opts.showDevParametrization  = 1;
  end
  
  
  if ~isfield(opts,'nMaxShow')
        opts.nMaxShow = .6;
  end
    
if  isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end
    
%% Make Fgiure

N = length(digdata.FileNames);
clear hFs
for kk=1:(ceil(N/pMax))    
    nStart=(kk-1)*pMax+1;
    nEnd=min([pMax*kk N]);
    
    hFs(kk)=figure('Name', ['Radial ' num2str(kk) ' '  FigLabel], 'Visible', 'On', ...
        'NumberTitle','off','color','w','units','pixels'); 
    hF=hFs(kk);
    hF.Position(1)=20;
    hF.Position(2)=50;
    hF.Position(3)=1850;
    hF.Position(4)=950;
    clf;

    t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left');
    t.Position(4)=t.Extent(4);
    t.Position(3)=hF.Position(3);
    t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

    for ii=nStart:nEnd
        % Create axes object
        ax=axes('parent',hF,'units','pixels');
        set(ax,'FontSize',8,'XMinorTick','on','YMinorTick','on',...
            'Box','on','YGrid','on','XGrid','on','units','pixels',...
            'yaxislocation','left');
        hold on;     
        cla;
        
        index = ii-nStart+1;
        n = nEnd-nStart+1;
        [a,b,c,d]=getAxesPos(index,n,hF.Position(3),hF.Position(4));
        ax.Position = [a b c d];   
        
        errorbar(digdata.r(:,ii),digdata.Zr(:,ii),...
            digdata.Zr_std(:,ii)./sqrt(digdata.nr(:,ii)),'linewidth',1);
        xlim([0 opts.rMaxShow]);    
        ylim([0 opts.nMaxShow]);
        xlabel('radial position (sites)');
        ylabel('charge');      
                
        % Draw the analysis string box
        
        % x variable string
        str = [digdata.FileNames{ii}];
        str = [str newline digdata.xVar ' : '];
        
        if isequal(digdata.xVar,'ExecutionDate')
            str = [str datestr(digdata.X(ii),'mm/DD HH:MM:ss')];
        else
            str = [str num2str(digdata.X(ii))];
        end        
        str = [str newline num2str(digdata.Natoms(ii)) ' atoms'];
        
        str = [str newline char(916) 'r:' num2str(digdata.rBinStep) 'sites'];

        % Draw the iteration number and variable value
        text(3, ax.Position(4)-2, ...
            str, ...
            'Units', 'pixels',...
            'FontSize', 8,...
            'verticalalignment','cap','HorizontalAlignment','left',...
            'interpreter','none');
        
        if opts.showDevParametrization
            w=ax.Position(3)*.3;
            h=ax.Position(4)*.3;
            axSub = axes('parent',hF,'units','pixels','fontsize',6);
            axSub.Position(1)  = [ax.Position(1)+ax.Position(3)-w];
            axSub.Position(2)  = [ax.Position(2)+ax.Position(4)-h];
            axSub.Position(3)  = w;
            axSub.Position(4)  = h;
            
            plot(digdata.Zr(:,ii),digdata.Zr_std(:,ii),'o');  
            tt=linspace(0,max(digdata.Zr(:,ii)),1e3);
            hold on
            plot(tt,sqrt(tt.*(1-tt)),'r-');
            
            text(.99,.01,'charge','units','normalized','fontsize',6,...
                'horizontalalignment','right',...
                'verticalalignment','bottom');
            text(.01,.99,['charge' newline 'sigma'],'units','normalized','fontsize',6,...
                'horizontalalignment','right',...
                'verticalalignment','top','Rotation',90);
        end        
    end      
    
    
    
    disp('done.');
    
end
end

function [axX,axY,axWidth,axHeight]=getAxesPos(nInd,nTot,xSize,ySize)
nInd=nInd-1;
yTop=30;
yBot=30;

xLeft=35;
xRight=20;

ySpace=35;
xSpace=35;

nRow=ceil(sqrt(nTot));

axHeight=(ySize-yTop-yBot-ySpace*(nRow-1))/nRow;
axWidth=(xSize-xLeft-xRight-xSpace*(nRow-1))/nRow;

axX=xLeft+(axWidth+xSpace)*mod(nInd,nRow);
axY=(ySize-yTop-axHeight)-floor(nInd/nRow)*(axHeight+ySpace);
end


