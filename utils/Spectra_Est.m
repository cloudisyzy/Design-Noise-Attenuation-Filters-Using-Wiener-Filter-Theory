function [pxx, f] = Spectra_Est(x, Option, plotFlag)
%   Using Welch Method to estimate the spectrum of given signals and give
%   plot
%   x - signal to be estimated
%   Option ='full', or 'half'
%   plotFlag =1, plot; =0, do not plot
%   pxx - estimated spectrum
%   f - nomarlized freq, from -0.5 to 0.5
%   Author: Ziyue Yang
%   Date: 2024.01.27

    [pxx, f] = pwelch(x, hamming(64), 32, 2048, 1, 'twosided');
    N = length(f);
    halfN = floor(N/2);

    if strcmp(Option, 'full')
        f = f - 0.5;
        pxx = [pxx(halfN+1:end); pxx(1:halfN)];
    elseif strcmp(Option, 'half')
        f = f(1:halfN);
        pxx = pxx(1:halfN);
    end

    if plotFlag
        semilogy(f, pxx, 'LineWidth',1.5)
    end

end

