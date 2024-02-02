function [p_AR, H_AR, f] = Spectra_AR(A, sigma2, Option, plotFlag)
%   Compute theoretical PSD of AR model given
%   A - denominator A(q)
%   sigma2 - variance of AR input noise
%   Option ='full', or 'half'
%   plotFlag - =1, plot; =0, do not plot
%   p_AR - theoretical PSD for given AR parameters
%   H_AR - frequency response of AR transfer func
%   f - nomarlized freq, from -0.5 to 0.5
%   Author: Ziyue Yang
%   Date: 2024.01.27

    f = linspace(-0.5, 0.5, 2048);
    H_AR = freqz(1, A, 2*pi*f);
    p_AR = abs(H_AR).^2 * sigma2;

    if strcmp(Option, 'full')

    elseif strcmp(Option, 'half')
        N = length(f);
        halfN = floor(N/2);
        f = f(halfN+1:end);
        H_AR = H_AR(halfN+1:end);
        p_AR = p_AR(halfN+1:end);
     
    end
    
    if plotFlag
        semilogy(f, p_AR, 'b--', 'LineWidth',1.5)   
    end

end

