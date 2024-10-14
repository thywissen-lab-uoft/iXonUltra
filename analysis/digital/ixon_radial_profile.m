function [Tics,Average,dev,n]=ixon_radial_profile(data,radial_step)
%main axii cpecified:
x=(1:size(data,2))-size(data,2)/2;
y=(1:size(data,1))-size(data,1)/2;
% coordinate grid:
[X,Y]=meshgrid(x,y);
% creating circular layers
Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % illustrating the principle:
% % figure;imagesc(Z_integer.*data)
% very fast MatLab calculations:
Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
Average=accumarray(Z_integer(:),data(:),[],@mean);

max_z = max(Z_integer(:));
figure(5);
Tics(1:5)
for kk=1:5
    zme = data(Z_integer(:)==kk);
    zme(zme<=100)=[];
    % keyboard
     [n,edges,bin] =histcounts(zme,30);
     centers = (edges(2)-edges(1))+edges;
     centers(end)=[];

     cmap = jet(5);
    
     plot(centers,n,'o','color',cmap(kk,:),'markerfacecolor',cmap(kk,:));
    hold on
end

dev=accumarray(Z_integer(:),data(:),[],@std);

n= accumarray(Z_integer(:),data(:),[],@(x) numel(x));

end