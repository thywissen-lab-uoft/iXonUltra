% function out = dig_Fidelity(Zdig,n1,n2,opts)
function digdata = dig_Fidelity(digdata)
% Author : CJ Fujiwara
%
% This codes calculates the fidelity of digitized data

% Should there be three figures?  
% (1) Single Shot FIdelity 
% (2) Average image Fidelity
% (3) Summary of Fidelity

if size(digdata.Zdig,4)==1
    warning('cannot compute fidelity on single shot image')
    return;
end

%% Calculate Fidelity
for kk = 1:length(digdata.FileNames)
    z1 = digdata.Zdig(:,:,kk,1);    % Image 1
    z2 = digdata.Zdig(:,:,kk,2);    % Image 2
    N1 = sum(z1,'all');             % Atom Number 1
    N2 = sum(z2,'all');             % Atom Number 2
    dz = z1 - z2;                   % Differential Image
    Nlost = N1 - N2;                % Number Lost
    Rlost = Nlost/N1;               % Percentage Loss
    Nhop = sum(abs(dz),'all')-Nlost;% Number Hopped
    Rhop = Nhop/N1;                 % Percentage hop


    digdata.lost_number(kk)       = Nlost;
    digdata.lost_fraction(kk)     = Rlost;
    digdata.hop_number(kk)        = Nhop;
    digdata.hop_fraction(kk)      = Rhop;
end

end
    
  
    
    
    % plot(xc,yc,'o','color','r','markersize',3,'markerfacecolor','r');
    
    % 
    % cc = parula(length(edges));
    % N = length(edges);
    % 
    % % rList = [20 40 60 80];
    % tt=linspace(0,2*pi,100);
    % for kk=2:length(edges)
    %     plot(xc+edges(kk)*cos(tt),yc+edges(kk)*sin(tt),'color',cc(kk,:))
    % end
    % xlim([min(n1) max(n1)]);
    % ylim([min(n2) max(n2)]);
    % 
    % % subplot(4,4,[11 12]);
    % subplot(2,4,[4]);
    % 
    % histogram(Revent,edges);
    % xlabel('radial position (sites)')
    % ylabel('loss or hop occurence')
    % xlim([0 max(edges)]);
    % set(gca,'box','on','linewidth',1,'fontname','times','fontsize',8);
    % 
    % % subplot(4,4,[15 16]);
    % subplot(2,4,[8]);
    % 
    % co=get(gca,'colororder');
    % xlabel('radial position (sites)')
    % % plot(r,N_expect,dev,'o','markerfacecolor',co(1,:),'linewidth',1,'markersize',6,...
    %     % 'color',co(1,:)*.5);
    %     % keyboard
    %     plot(r(1:N),N_expect(1:N),'k-');
    %     hold on
    % 
    % ps=scatter(r(1:N),N_expect(1:N),40,cc,'filled');
    % set(ps,'markeredgecolor','k')
    % ylabel('average occupation')
    % xlim([0 max(edges)]);
    % xlabel('radial position (sites)')
    % set(gca,'box','on','linewidth',1,'fontname','times','fontsize',8);
    %%

% end


% end

% function [Tics,Average,dev,n]=radial_profile(data,radial_step)
% %main axii cpecified:
% x=(1:size(data,2))-size(data,2)/2;
% y=(1:size(data,1))-size(data,1)/2;
% % coordinate grid:
% [X,Y]=meshgrid(x,y);
% % creating circular layers
% Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % % illustrating the principle:
% % % figure;imagesc(Z_integer.*data)
% % very fast MatLab calculations:
% Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
% Average=accumarray(Z_integer(:),data(:),[],@mean);
% 
% 
% dev=accumarray(Z_integer(:),data(:),[],@std);
% 
% n= accumarray(Z_integer(:),data(:),[],@(x) numel(x));
% 
% end