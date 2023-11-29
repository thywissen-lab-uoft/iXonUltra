function data = ixon_PositionHistogram(data)
fprintf('Calculating position space histogram')

for n=1:length(data)
    ROI = data(n).ROI;
    c1 = find(data(n).X>=ROI(1),1);
    c2 = find(data(n).X>=ROI(2),1);
    r1 = find(data(n).Y>=ROI(3),1);
    r2 = find(data(n).Y>=ROI(4),1);
    x = c1:c2;
    y = r1:r2; 

    z = data(n).ZNoFilter(y,x,1);

    if isfield(data(n),'RotationMask')
        m = data(n).RotationMask;
        z = z(m);
    end

    % Perform Histogram on first no filter image to get the edges
    [~,xe] = histcounts(z,1000);

    for kk=1:size(data(n).ZNoFilter,3)          
        z = data(n).ZNoFilter(y,x,kk);
        if isfield(data(n),'RotationMask')
            z = z(m);
        end        
        [N,edges] = histcounts(z,xe);
        centers = (edges(1:end-1) + edges(2:end))/2;
        data(n).HistogramNoFilter(kk) = ...
            struct('Edges',edges,'Centers',centers,'N',N);
    end   

    for kk=1:size(data(n).Z,3)    
        z = data(n).Z(y,x,kk);
        if isfield(data(n),'RotationMask')
            z = z(m);
        end 
        [N,edges] = histcounts(z,xe);
        centers = (edges(1:end-1) + edges(2:end))/2;
        data(n).Histogram(kk) = ...
            struct('Edges',edges,'Centers',centers,'N',N);
    end 
end

disp('done')


end

