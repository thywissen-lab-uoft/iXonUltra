function [hF] = dig_showCloud(digdata,opts)

if nargin ==1
    opts = struct;
end

if ~isfield(opts,'doAnimate')
    opts.doAnimate = 0;
end

if ~isfield(opts,'doSave')
    opts.doSave = 0;
end

if ~isfield(opts,'ROI')
    opts.ROI = 'max';
end

%% Get ROI
    
    n1 = digdata.n1;
    n2 = digdata.n2;
    
    if isequal(opts.ROI,'max')
       opts.ROI = [min(n1) max(n1) min(n2) max(n2)]; 
    end    
    R = opts.ROI;
   
    
    %% Initialize Graphics
    hF = figure;
    clf
    hF.Color='w';
    hF.Position= [100 100 500 400];
    hF.Name = 'Digital';
    
    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        tFig=uicontrol('style','text','string',opts.FigLabel,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left');
        tFig.Position(4)=tFig.Extent(4);
        tFig.Position(3)=hF.Position(3);
        tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
    end    
    
  
 

    % Image
    ax=axes;
    hImg = imagesc(n1,n2,mean(digdata.Zdig,3));    
    rectangle('Position',[R(1) R(3) R(2)-R(1) R(4)-R(3)],...
        'EdgeColor','r')
    xlabel('n1 sites');
    ylabel('n2 sites');
    c=colorbar;
    c.Label.String = 'average counts/site';
    axis equal tight
    set(gca,'box','on','linewidth',1,'fontsize',10,'ydir','normal');
    ca = [0 0 0];       
    cb = [0.7 .1 .6];
    cc = [linspace(ca(1),cb(1),1000)' ...
        linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
    colormap(hF,cc);
    colormap(hF,'bone');

    c.Label.String = 'counts/site';
    caxis([0 1])

       % String Labels
    tNW = text(0.01,.98,'','units','normalized','verticalalignment',...
        'top','horizontalalignment','left','fontsize',8,'Color','r');

    tSW = text(0.01,0.01,'','units','normalized','verticalalignment',...
        'bottom','horizontalalignment','left','fontsize',8,'color','r',...
        'Visible','off','interpreter','none');

    tSE = text(.99,.01,'','units','normalized','verticalalignment',...
        'bottom','horizontalalignment','right','fontsize',8,'color','r',...
        'Visible','off');
    %% Populate with data
    
    if opts.doAnimate
        for nn = 1:size(digdata.Zdig,3)

            set(hImg,'Cdata',digdata.Zdig(:,:,nn));

            if isequal(digdata.xVar,'ExecutionDate')
                sX = datestr(digdata.X(nn),'yy-mm-dd_HH-MM-SS');
            else
                sX = num2str(digdata.X(nn));
            end
            sSW = [digdata.xVar newline sX];
            set(tSW,'String',sSW,'Visible','on');

            sNW = [num2str(nn) '/' num2str(size(digdata.Zdig,3))];
            set(tNW,'String',sNW,'Visible','on');

            sSE = ['threshold : ' num2str(digdata.Threshold(nn)) newline ...
                num2str(sum(digdata.Zdig(:,:,nn),'all')) ' atoms'];
            set(tSE,'String',sSE,'Visible','on');   
            
            frame=getframe(hF);
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);              
            filename = fullfile(opts.saveDir,opts.filename);            
            switch nn
                case 1
                    imwrite(A,map,filename,'gif','LoopCount',...
                        Inf,'DelayTime',1);
                case length(digdata)
                    imwrite(A,map,filename,'gif','WriteMode',...
                        'append','DelayTime',1);
                otherwise
                    imwrite(A,map,filename,'gif','WriteMode',...
                        'append','DelayTime',.1);
             end
               
        
            
        end

        
    end
    
    %% Last Image
    
              
            c.Label.String = 'average counts/site';

            set(hImg,'Cdata',mean(digdata.Zdig,3));
%                         set(hImg,'Cdata',imgaussfilt(mean(digdata.Zdig,3),.5));

            caxis(ax,'auto')
            % tSE.Visible = 'off';
            tSW.Visible = 'off';

            set(tNW,'String',[num2str(size(digdata.Zdig,3)) ' images']);

            sSE = ['avg threshold :' num2str(mean(digdata.Threshold)) newline ...
                num2str(round(mean(sum(digdata.Zdig,[1 2])))) '\pm' ...
                num2str(round(std(sum(digdata.Zdig,[1 2])))) ' atoms'];

            set(tSE,'String',sSE);

end













