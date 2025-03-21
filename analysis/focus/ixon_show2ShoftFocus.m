function [hF] = ixon_show2ShoftFocus(focus,opts)

if nargin ==1
    opts=struct;
end

for kk=1:numel([focus.ExecutionDate])

    if ~isfield(opts,'Parent') || isempty(opts.Parent)
        hF(kk) = figure;
        set(hF(kk),'color','w');
        hF(kk).Position=[50 50 600 550];
    else
        hF = opts.Parent;
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
        ax(jj) = subplot(m,n,jj,'parent',hF(kk));
        imagesc(focus.Images_Pos(:,:,kk,jj)*focus.Counts(kk,jj),'parent',ax(jj));
        axis equal tight
        set(ax(jj),'XTickLabel',{},'YTickLabel',{});
        title(ax(jj),['Image ' num2str(jj) ' : ' ...
            num2str(focus.Piezos(kk,jj)) ' V']);

        axFreq = subplot(m,n,4,'parent',hF(kk));
        ps(jj) = plot(focus.RadialFrequency,...
            focus.RadialFFT(:,kk,jj),'.-',...
            'parent',axFreq);  
        hold on
        s=['(' num2str(jj) ') ' num2str(focus.Piezos(kk,jj)) ' V, ' ...
            num2str(focus.Scores(kk,jj),'%.1f') '/px'];
        legStr{jj}=s;
        c=get(ax(jj),'CLim');
        cL(jj)=c(2);
        text(.01,.01,[num2str(focus.Counts(kk,jj),'%.2e') ' counts'],...
            'units','normalized','verticalalignment','bottom',...
            'parent',ax(jj),'color','w');
    end
    
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

