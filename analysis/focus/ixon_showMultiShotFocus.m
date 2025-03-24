function [hF] = ixon_showMultiShotFocus(focus,opts)

if nargin ==1
    opts=struct;
end

for kk=1:size([focus.Scores],1)


  if ~isfield(opts,'Parent')
    FigName = 'Focus';      
    fig=figure;
    fig.Color='w';
    fig.Name=FigName;
    % fig.ToolBar='none';
    % fig.MenuBar='none';
    fig.Position=[5 50 600 500];
    clf(fig);
    my_parent=fig;
    hF(kk)=my_parent;

    else
        fig = opts.Parent;
        for ll=1:length(fig.Children)
            delete(fig.Children(1))
        end
        my_parent=fig;
        hF=my_parent;
    end



    N = size(focus.Images_Pos,4);
    legStr={};
    cL=[];
    if size(focus.Scores,2)==3
        m = 2;
        n = 2;
    else
        m = N+1;
        n = 1;
    end

    for jj=1:N
        ax(jj) = subplot(m,n,jj,'parent',my_parent);   
        
        imagesc(focus.Images_Pos(:,:,kk,jj)*focus.Counts(kk,jj),'parent',ax(jj));
        axis(ax(jj),'equal');
        axis(ax(jj),'tight');
        set(ax(jj),'XTickLabel',{},'YTickLabel',{});
        title(ax(jj),['Image ' num2str(jj) ' : ' ...
            num2str(focus.Piezos(kk,jj)) ' V']);

        axFreq = subplot(m,n,4,'parent',my_parent);
        ps(jj) = plot(focus.RadialFrequency,...
            focus.RadialFFT(:,kk,jj),'.-',...
            'parent',axFreq);  
        hold(axFreq,'on');
        s=['(' num2str(jj) ') ' num2str(focus.Piezos(kk,jj)) ' V, ' ...
            num2str(focus.Scores(kk,jj),'%.1f') '/px'];
        legStr{jj}=s;
        c=get(ax(jj),'CLim');
        cL(jj)=c(2);
        text(.01,.01,[num2str(focus.Counts(kk,jj),'%.2e') ' counts'],...
            'units','normalized','verticalalignment','bottom',...
            'parent',ax(jj),'color','k');
    end

    uicontrol('parent',my_parent,'style','text','backgroundcolor','w',...
        'Position',[1 1 300 20],'string',focus.Params(kk).ExecutionDateStr,...
        'horizontalalignment','left')
    
    for jj=1:N
        set(ax(jj),'CLim',[0 max(cL)]);
    end
    voff=focus.Params(kk).piezo_offset;
    text(.99,.01,[num2str(voff) 'V offset'],'parent',axFreq,...
            'units','normalized','horizontalalignment','right',...
            'verticalalignment','bottom');
    str=['$F(k):=\int\mathrm{d}\phi f(k,\phi)$'];
    xlim(axFreq,[0 .5]);
    ylabel(axFreq,str,'interpreter','latex')
    title(axFreq,'Radial FFT')
    legend(ps,legStr,'location','northeast')
end

end

