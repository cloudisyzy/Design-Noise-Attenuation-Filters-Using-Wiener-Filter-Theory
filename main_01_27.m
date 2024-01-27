clear
addpath('mfiles\')
addpath('skel\')
addpath('utils\')

% y <=> noisy speech; x <=> clean speech; v <=> additive noise
% y = x + v
% z = [y,v,y,v,...]
[z,fs] = audioread('EQ2401project1data2024.wav');
% Extract noisy speech samples 
y = z([4860:18780 26660:41212 51340:end]); % y = z([4501:18999 26001:41999 51001:end]);
% soundsc(y, fs)
% Extract pure noise samples 
v = z([1:4700, 18900:26500, 41450:51000]); % v = z([1:4500, 19000:26000, 42000:51000]);
% soundsc(v, fs)
% Estimated speech samples 
x = y - [v; v; v(1:length(y)-2 * length(v))];

N_v = 11; % 19 by visual inspection; 11 by computation
N_y = 17; % 30 by visual inspection; 17 by computation
N_x = N_y;

%% Compute auto-correlations

L = 100;
r_yy = xcovhat(y, y, L);
r_vv = xcovhat(v, v, L);
r_xx = r_yy - r_vv;
r_yx = r_xx;
R_yy = covhat(y, L);

%% Spectrum Estimation

[A_v, sigma2_v] = levinson(r_vv, N_v);
% [A_v, sigma2_v] = ar_id(v, N_v);
[A_x, sigma2_x] = levinson(r_xx, N_x);
% [A_x, sigma2_x] = ar_id(x, N_x);
[A_y, sigma2_y] = levinson(r_yy, N_y);
% [A_y, sigma2_y] = ar_id(y, N_y);

[PhixyNum,PhixyDen,PhiyyNum,PhiyyDen] = spec_add(A_x, sigma2_x, A_v, sigma2_v);

%% FIR

n_fir = 30;
SigmaYy = r_yy(1:n_fir);
SigmaVv = r_vv(1:n_fir);
SigmaXx = SigmaYy - SigmaVv;
SigmaYx = SigmaXx;
SigmaYY = R_yy(1:n_fir, 1:n_fir);

[xhatfir, thetahatfir] = firw(z, SigmaYx, SigmaYY);

%% Non-Causal

[xhatnc, numnc, dennc] = ncw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen);

%% Causal

% m=0, filtering
m = 0;
[xhatc, numc, denc] = cw(z, PhixyNum, PhixyDen, PhiyyNum, PhiyyDen, m);

%% Plots

figure;
spec_comp(A_x, sigma2_x, A_v, sigma2_v, numnc, dennc, numc, denc, thetahatfir);

figure;
t_c = [0:30]';
s_c = (t_c==0);
h_c = filter(numc, denc, s_c);
t_nc = [-30:30]';
s_nc = (t_nc==0);
h_nc = ncfilt(numnc, dennc, s_nc);
stem(0:length(thetahatfir)-1, thetahatfir, "+", LineWidth=2)
hold on
stem(t_c', h_c, "x", LineWidth=2)
hold on
stem(t_nc', h_nc, "o", LineWidth=2)
title('Impulse Response of Three Filters')
xlabel('Time Lag: n')
ylabel('h[n]')
legend('FIR', 'Causal', 'Non-Causal')

figure;
[~,~,~] = Spectra_AR(A_v, sigma2_v, 'half', 1);hold on
[~,~] = Spectra_Est(v, 'half', 1);
title('Noise Order Estimation via Plot')
xlabel('Normalized Frequency \nu, unit:Hz')
ylabel('Power Spectrum, unit:dB')
legend(sprintf('PSD of AR-%d', N_v), 'Estimated Spectrum of Noise'); hold off

figure;
[~,~,~] = Spectra_AR(A_y, sigma2_y, 'half', 1); hold on
[~,~] = Spectra_Est(y, 'half', 1);
title('Noisy Speech Order Estimation via Plot')
xlabel('Normalized Frequency \nu, unit:Hz')
ylabel('Power Spectrum, unit:dB')
legend(sprintf('PSD of AR-%d', N_y), 'Estimated Spectrum of Noisy Speech'); hold off

% % the correctness of the below estimation need to be testified
% figure;
% [~,~,~] = Spectra_AR(A_x, sigma2_x, 'half', 1); hold on
% [~,~] = Spectra_Est(x, 'half', 1);
% title('Speech Order Estimation via Plot')
% xlabel('Normalized Frequency \nu, unit:Hz')
% ylabel('Power Spectrum, unit:dB')
% legend(sprintf('PSD of AR-%d', N_x), 'Estimated Spectrum of Speech'); hold off

%% Play Sound

delta_t = length(z)/fs + 0.5;
soundsc(z, fs)
pause(delta_t)
soundsc(xhatfir, fs)
pause(delta_t)
soundsc(xhatnc, fs)
pause(delta_t)
soundsc(xhatc, fs)

%% Save Audio Files

% audiowrite('results\fir_visual.wav', xhatfir, fs);
% audiowrite('results\nc_visual.wav', xhatnc, fs);
% audiowrite('results\c_visual.wav', xhatc, fs);
% audiowrite('results\_origin.wav', z, fs);

% pause_data = zeros(10000, 1);
% blind_test = [xhatc; pause_data; z; pause_data; xhatfir; pause_data; xhatnc]; % mse
% blind_test = [xhatnc; pause_data; xhatfir; pause_data; z; pause_data; xhatc]; % visual
% audiowrite('results\blind_test_visual.wav', blind_test, fs);