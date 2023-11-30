function [output] = bimodalPDFFit(z)
    z=z(:);
    z(isnan(z))=[];
    % Divide into 2 clusters (n=0 vs. n=1)
    [idx,c]=kmeans(z,2);
    z1 = z(idx==1); % No atom cluster
    z2 = z(idx==2); % atom cluster

    % The n=0 (no atom) distribution is fit with a gaussian + a gamma
    % function.  The physical motivation is that the gaussian represents
    % the actual noise of measurement (read + subtraction) while the gamma
    % tail represents counts from sites filled with atoms leaking into
    % unoccupied sites.
    try 
        % Find numerical PDF for n=0 atoms;
        [f,zlist] = ksdensity(z1);
        [~,ip] = max(f);
        zp = zlist(ip);

        % Divide the n = 0 into to halves
        z1_left = z1(z1<zp);
        z1_right = z1(z1>=zp);

        % Fit the left side of the distribution to a half gaussian
        z1_left_flipped = max(z1_left)-z1_left;
        pd1_left = fitdist(z1_left_flipped,'Half Normal');
        s1_left = pd1_left.sigma;
        mu1_left = zp;
        

        % Fit the Right side to a gamma function
        pdg=fitdist(z1_right,'gamma');

        % Create Combined PDF
        warning off
        LB = [mu1_left-s1_left s1_left/2 0 0 .3];
        [pdf0_c,pdf0_cint] = mle(z1,'pdf',@pdf_gauss_gamma,...
            'start',[mu1_left s1_left pdg.a pdg.b 0.5],'alpha',.1,...
            'LowerBound',LB);
        warning on
        
    catch ME

    end
        
        
    % The n=1 (atom) distribution is fit with a simple gaussian
    % distribution
    % pd1 = fitdist(z2,'normal'); % This works too, but use MLE to make it
    % the same output
    [pdf1_c,pdf1_cint] = mle(z2,'distribution','normal');

    % Create output
    output = struct;

    % n = 0 probability distribution function
    output.pdf0 = @(x) pdf_gauss_gamma(x,pdf0_c(1),pdf0_c(2),pdf0_c(3),pdf0_c(4),pdf0_c(5));
    output.pdf0_coeffs = pdf0_c;
    output.pdf0_cints  = pdf0_cint;
    output.pdf0_func = @pdf_gauss_gamma;
    output.pdf0_counts = z1;
    output.pdf0_centroid = c(1);


    % n = 1 probability distribution function
    output.pdf1 = @(x) normpdf(x,pdf1_c(1),pdf1_c(2));
    output.pdf1_coeffs = pdf1_c;
    output.pdf1_cints  = pdf1_cint;
    output.pdf1_func = @normpdf;
    output.pdf1_counts = z2;
    output.pdf1_centroid = c(2);


end

function y = pdf_gauss_gamma(x,u1,s1,a,b,alpha)
    y = (alpha*gampdf(x,a,b)+(1-alpha)*normpdf(x,u1,s1));
end

