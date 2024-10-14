function bin_showRadialHistogram(bindata,opts)

if nargin ==1
    opts= struct;
end


if ~isfield(opts,'FigLabel')
    opts.FigLabel=[];
end
    if ~isfield(opts,'doAnimate')
       opts.doAnimate = 0; 
    end
%% Initialize Graphics
if isfield(opts,'FigureNumber')
    hF = figure(opts.FigureNumber);
    clf
else    
    hF = figure;
end

    hF.Color='w';
    hF.Position= [0 50 800 800];
    hF.Name = 'Binned Histogram';
    
    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        tFig=uicontrol('style','text','string',opts.FigLabel,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left');
        tFig.Position(4)=tFig.Extent(4);
        tFig.Position(3)=hF.Position(3);
        tFig.Position(1:2)=[5 hF.Position(4)-tFig.Position(4)];
    end    

    markers={'o','s','t','d'};

    for ii=1:length(bindata)
        clf
        J = length(bindata(ii).LatticeRadialHistogram);
        % J = 1;
        pfits=[];
        fit_strs={};
        for jj=1:J

            r = bindata(ii).LatticeRadialHistogram(jj).RadialVector;
            c = bindata(ii).LatticeRadialHistogram(jj).Centers;
        
            nmax = 6;
            cmap = jet(nmax);

            % nMax = max([])
            % for rr = 1:size(bindata(ii).LatticeRadialHistogram(jj).N,1)
            for rr=1:nmax
                
                ax = subplot(nmax,2,2*rr-1);

                N = bindata(ii).LatticeRadialHistogram(jj).N(rr,:);
                plot(c,N,[markers{jj}],'color',.4*cmap(rr,:),'linewidth',1,...
                    'markerfacecolor',cmap(rr,:),'markersize',8);
                hold on
                xlabel('counts');
                ylabel('occurences');
                xlim([1000 8000]);

                Nt = bindata(ii).LatticeRadialHistogram(jj).ClusterThreshold(rr);

                Nscaled = N/Nt;
                
            end

                subplot(2,2,2);
                r = bindata(ii).LatticeRadialHistogram(jj).RadialVector;
                t = bindata(ii).LatticeRadialHistogram(jj).ClusterThreshold;
    
                 [cface,cedge] = ixoncolororder(jj);
                plot(r,t,['-' markers{jj}],'color',cedge,'linewidth',2,...
                    'markerfacecolor',cface,'markersize',10);
                hold on
                set(gca,'yaxislocation','right');
                xlabel('radial position (sites)');
                ylabel('threshold (counts)');

                title('threshold')


                  if isfield(bindata(ii).LatticeRadialHistogram(jj),'GaussFit')
                    rrr=linspace(0,max(r),100);
                    pfits(end+1)=plot(rrr,feval(bindata(ii).LatticeRadialHistogram(jj).GaussFit,rrr),...
                        '-','linewidth',2,'color',cedge*.5);
                    fit_strs{end+1}=['$A:' num2str(round(bindata(ii).LatticeRadialHistogram(jj).Amplitude)) ...
                        ',~s:' num2str(round(bindata(ii).LatticeRadialHistogram(jj).Radius)) '$'];

                  end                

                subplot(2,2,4);
                 [cface,cedge] = ixoncolororder(jj);
                % plot(r,t,['-' markers{jj}],'color',cedge,'linewidth',2,...
                    % 'markerfacecolor',cface,'markersize',10);

                Zbsc = bindata(ii).LatticeRadialHistogram(jj).ZbinScale;
                Zbsc(Zbsc<=0)=[];
                [n,edges]=histcounts(Zbsc,30);
                c = edges(2)-edges(1)+edges;
                c(end)=[];

                plot(c,n,'o-','color',cedge,'linewidth',2,...
                    'markerfacecolor',cface,'markersize',10)

                hold on
                set(gca,'yaxislocation','right');
                xlabel('$\mathrm{counts}/N_{\mathrm{thresh}}(r)$',...
                    'interpreter','latex');
                ylabel('occurences');   
                % xlim([.75 3]);
                xlim([1000 8000])
                title('scaled histgoram')

        end

          if isfield(bindata(ii).LatticeRadialHistogram(1),'GaussFit')
                subplot(2,2,2);
                legend(pfits,fit_strs,'interpreter','latex');
          end

            
            

            
            for rr = 1:nmax
                subplot(nmax,2,2*rr-1);
                rMe=[];
                Nt=[];
                for jj=1:J
                    rMe(jj) = bindata(ii).LatticeRadialHistogram(jj).RadialVector(rr);
                    Nt(jj) = bindata(ii).LatticeRadialHistogram(jj).ClusterThreshold(rr);
                end

                rMe = round(rMe);
                Nt = round(Nt);

                srMe = sprintf('%.0f,' , rMe);
                srMe = srMe(1:end-1);% strip final comma


                sNt = sprintf('%.0f,' , Nt);
                sNt = sNt(1:end-1);% strip final comma

                str=['$r:[' srMe  ']$' newline ...
                    '$\mathrm{thresh}:[' sNt ']$'];
                text(.99,.99,str,'units','normalized','interpreter','latex',...
                    'verticalalignment','top','horizontalalignment','right',...
                    'fontsize',10);
            end



     
    end

end

