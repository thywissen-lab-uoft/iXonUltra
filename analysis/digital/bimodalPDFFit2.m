function [output] = bimodalPDFFit2(z)

    if sum(z<0,'all')>0
        warning on
        warning('negative counts detected. why is this happening?')
        z(z<0)=0;
    end

    z=z(:);
    z(isnan(z))=[];

    % Divide into 2 clusters (n=0 vs. n=1)
    [idx,c,sumD,D]=kmeans(z,2);
    
    z1 = z(idx==1); % No atom cluster
    z2 = z(idx==2); % atom cluster

    % The n=0 (no atom) distribution is fit with a gaussian + a gamma
    % function.  The physical motivation is that the gaussian represents
    % the actual noise of measurement (read + subtraction) while the gamma
    % tail represents counts from sites filled with atoms leaking into
    % unoccupied sites.
    
    

    % Fit the Right side to a gamma function
    [pdf0_c,pdf0_cint]=mle(z1,'distribution','Exponential');
        
    % The n=1 (atom) distribution is fit with a simple gaussian
    % distribution
    % pd1 = fitdist(z2,'normal'); % This works too, but use MLE to make it
    % the same output
    [pdf1_c,pdf1_cint] = mle(z2,'distribution','normal');

    
    % Calculate the fidelity via numerical integration
    try
        s1 = integral(@(x) exppdf(x,pdf0_c(1)),0,2*max(z));
        s2 = integral(@(x) normpdf(x,pdf1_c(1),pdf1_c(2)),0,2*max(z));
        f = 0.5*integral(@(x) abs(exppdf(x,pdf0_c(1))-normpdf(x,pdf1_c(1),pdf1_c(2))),0,2*max(z));       
        
    catch ME
        f = NaN;
        
    end

    % Create output
    output = struct;
    
    output.Fidelity = f;

    % n = 0 probability distribution function
    output.pdf0 = @(x) exppdf(x,pdf0_c(1));
    output.pdf0_coeffs = pdf0_c;
    output.pdf0_cints  = pdf0_cint;
    output.pdf0_func = @exppdf;
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

