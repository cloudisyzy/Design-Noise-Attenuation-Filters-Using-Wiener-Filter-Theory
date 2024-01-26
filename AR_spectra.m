clear

addpath("mfiles\")
addpath("skel\")

% Params
sigma2 = 1;
K = 10000;
A = [1,0.2,-0.4,-0.8];
e = sqrt(sigma2)*randn(K,1);
y = filter(1,A,e);

% Welch Estimated Spectrum
[pxx, f_welch] = pwelch(y, hamming(256), 128, 512, 1, 'twosided');
f_welch = f_welch - 0.5;
N = length(f_welch);
halfN = floor(N/2);
pxx = [pxx(halfN+1:end); pxx(1:halfN)];

figure;
plot(f_welch, 10*log10(pxx))
hold on
xlabel('Frequency (Normalized)')
ylabel('Spectral Density (Estimated)')
title(sprintf('Estimated vs True Spectrum of AR-%d Model', length(A)-1))

% Theoretical Spectrum
f_theory = linspace(-0.5, 0.5, 1000);
H = freqz(1, A, 2*pi*f_theory);
S_y = abs(H).^2 * sigma2;
plot(f_theory, 10*log10(S_y));
legend('Welch Estimation','Ground Truth')
hold off

% % Alternative: from 0 to 0.5 freq, must be correct
% sigma2 = 1;
% K = 10000;
% A = [1,0.2,-0.4,-0.8];
% e = sqrt(sigma2)*randn(K,1);
% y = filter(1,A,e);
% 
% figure;
% [pxx, f] = pwelch(y, hamming(256), 128, 512, 1);
% plot(f, 10*log10(pxx))
% xlabel('Frequency (Normalized)')
% ylabel('Spectral Density (Estimated)')
% title('Welch PSD Estimate of AR(1) Model')
% 
% figure;
% f = linspace(0, 0.5, 500);
% H = freqz(1, A, 2*pi*f);
% S_y = abs(H).^2 * sigma2;
% plot(f, 10*log10(S_y));
% xlabel('Frequency (Normalized)');
% ylabel('Spectral Density (True)');
% title('Theoretical Spectrum of AR(1) Model');

