clear 

addpath('mfiles\')
addpath('skel\')
[z,fs] = audioread('EQ2401project1data2024.wav');

% figure;
% plot(0:length(z)-1, z)

% y <=> noisy speech; x <=> clean speech; v <=> additive noise
% y = x + v

% Extract noisy speech samples 
y_1 = z(4860:18780);
y_2 = z(26660:41212);
y_3 = z(51340:68663);
y = [y_1; y_2; y_3];


% Extract pure noise samples 
v1 = z(1:4700);
v2 = z(18900:26500);
v3 = z(41450:51000);
v = [v1;v2;v3];
% sound(v, fs)

%% Estimate noise model (AR)
estimated_order_noise = 19;
% Plot the theoretical spectra and estimated noise spectra
[pxx, f_welch] = pwelch(v, hamming(64), 32, 1024, 1, 'twosided');
f_welch = f_welch - 0.5;
N = length(f_welch);
halfN = floor(N/2);
pxx = [pxx(halfN+1:end); pxx(1:halfN)];

[Anoisehat, sigma2noisehat] = ar_id(v, estimated_order_noise);
f_theory = linspace(-0.5, 0.5, 1000);
H = freqz(1, Anoisehat, 2*pi*f_theory);
S_y = abs(H).^2 * sigma2noisehat;

figure;
plot(f_welch, 10*log10(pxx))
hold on
plot(f_theory, 10*log10(S_y))
xlabel('Frequency (Normalized)')
ylabel('Spectral Density: dB')
title(sprintf('Noise Spectrum vs Spectrum of AR-%d Model', estimated_order_noise))
legend('Noise Spectra',sprintf('AR-%d Spectra', estimated_order_noise))

%% Estimate sound model (AR)

estimated_order_sound = 17;

e = sqrt(sigma2noisehat) * randn(length(y), 1);
noise_in_speech = filter(1, Anoisehat, e);
x_denoise = y - noise_in_speech;

[pxx, f_welch] = pwelch(x_denoise, hamming(64), 32, 1024, 1, 'twosided');
f_welch = f_welch - 0.5;
N = length(f_welch);
halfN = floor(N/2);
pxx = [pxx(halfN+1:end); pxx(1:halfN)];

[Ahat, sigma2hat] = ar_id(x_denoise, estimated_order_sound);
f_theory = linspace(-0.5, 0.5, 1000);
H = freqz(1, Ahat, 2*pi*f_theory);
S_y = abs(H).^2 * sigma2hat;

figure;
plot(f_welch, 10*log10(pxx))
hold on
plot(f_theory, 10*log10(S_y))
xlabel('Frequency (Normalized)')
ylabel('Spectral Density: dB')
title(sprintf('Speech Spectrum vs Spectrum of AR-%d Model', estimated_order_sound))
legend('Speech Spectra',sprintf('AR-%d Spectra', estimated_order_sound))

%% FIR
n_fir = 30; % FIR order
SigmaYyhat = xcovhat(y, y, n_fir);
SigmaVvhat = xcovhat(v, v, n_fir);
% R_yx = R_xx = R_yy - R_vv, since independency of noise
SigmaYxhat = SigmaYyhat - SigmaVvhat;
SigmaYYhat = covhat(y, n_fir);
[xhatfir, thetahatfir] = firw(z, SigmaYxhat, SigmaYYhat);

% figure;
% plot(xhatfir)
% hold on
% plot(z)
% legend('FIR', 'Noisy')

%% Non-Causal

[PhixyNum,PhixyDen,PhiyyNum,PhiyyDen] = ...
               spec_add(Ahat, sigma2hat, Anoisehat, sigma2noisehat);
[xhatnc, numnc, dennc] = ncw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen);
soundsc(xhatnc, fs)

%% Play Sound
audioLength = length(z) / fs;
% soundsc(z, fs);
pause(audioLength + 0.5);
% soundsc(xhatfir, fs)

[compare,fs] = audioread('赋值nc()30,0.041,0.025.wav');
soundsc(compare, fs)
% 
% pause(audioLength + 0.5);
% soundsc(z, fs);

mae(xhatnc,z)

