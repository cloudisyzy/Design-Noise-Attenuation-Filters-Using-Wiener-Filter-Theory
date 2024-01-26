clear 
[z,fs] = audioread('EQ2401project1data2024.wav');

estimated_order = 19;
% sound(z,fs)
% figure;
% plot(z)

z1 = z(1:4700);
z2 = z(18900:26500);
z3 = z(41450:51000);
z = [z1;z2;z3];

% [pxx, f_welch] = pwelch(z, hamming(64), 32, 512, 1, 'twosided');
[pxx, f_welch] = pwelch(z, hamming(64), 32, 1024, 1, 'twosided');
f_welch = f_welch - 0.5;
N = length(f_welch);
halfN = floor(N/2);
pxx = [pxx(halfN+1:end); pxx(1:halfN)];

[Ahat, sigma2hat] = ar_id(z, estimated_order);
f_theory = linspace(-0.5, 0.5, 1000);
H = freqz(1, Ahat, 2*pi*f_theory);
S_y = abs(H).^2 * sigma2hat;

figure;
plot(f_welch, 10*log10(pxx))
hold on
plot(f_theory, 10*log10(S_y))
xlabel('Frequency (Normalized)')
ylabel('Spectral Density: dB')
title(sprintf('Noise Spectrum vs Spectrum of AR-%d Model', length(Ahat)-1))
legend('Noise Spectra',sprintf('AR-%d Spectra', length(Ahat)-1))

