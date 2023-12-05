function pdfValue = skewNormalPDF(x, theta, alpha, omega) 
    % Skew Normal PDF implementation 
    % x: random variable 
    % theta: location parameter 
    % alpha: shape parameter controlling skewness 
    % omega: scale parameter 
 
    % Standard normal PDF and CDF 
    phi = @(x) exp(-0.5 * x.^2) / sqrt(2 * pi); 
    Phi = @(x) 0.5 * (1 + erf(x / sqrt(2))); 
 
    % Calculate the skew normal PDF 
    z = (x - theta) / omega;  % Standardize the random variable 
    pdfValue = 2 / omega * phi(z) .* Phi(alpha * z); 
end 