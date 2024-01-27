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

% Compute auto-correlations
L = 100;
r_yy = xcovhat(y, y, L);
r_vv = xcovhat(v, v, L);
r_xx = r_yy - r_vv;
r_yx = r_xx;
R_yy = covhat(y, L);

%% Compute the Best AR Order that Fits Noise/Speech

Nmax = 40; % Maximum test order

err_v_est = zeros(1, 40);
err_y_est = zeros(1, 40);

for model_order = 1:Nmax
    [A_v_temp, sigma2_v_temp] = levinson(r_vv, model_order);
    [p_v_AR,~,~] = Spectra_AR(A_v_temp, sigma2_v_temp, 'half', 0);
    [p_v,~] = Spectra_Est(v, 'half', 0);
    err_v_est(model_order) = mse(p_v, p_v_AR');

    [A_y_temp, sigma2_y_temp] = levinson(r_yy, model_order);
    [p_y_AR,~,~] = Spectra_AR(A_y_temp, sigma2_y_temp, 'half', 0);
    [p_y,~] = Spectra_Est(y, 'half', 0);
    err_y_est(model_order) = mse(p_y, p_y_AR');
end

[MMSE_v, MMSE_order_v] = min(err_v_est);
[MMSE_y, MMSE_order_y] = min(err_y_est);

%% Plot

figure;
[A_v, sigma2_v] = levinson(r_vv, MMSE_order_v);
[~,~,~] = Spectra_AR(A_v, sigma2_v, 'half', 1); hold on
[~,~] = Spectra_Est(v, 'half', 1);
title('Noise Order Estimation via Plot')
xlabel('Normalized Frequency \nu, unit:Hz')
ylabel('Power Spectrum, unit:dB')
legend(sprintf('PSD of AR-%d', MMSE_order_v), 'Estimated Spectrum of Noise');

figure;
[A_y, sigma2_y] = levinson(r_yy, MMSE_order_y);
[~,~,~] = Spectra_AR(A_y, sigma2_y, 'half', 1); hold on
[~,~] = Spectra_Est(y, 'half', 1);
title('Noisy Speech Order Estimation via Plot')
xlabel('Normalized Frequency \nu, unit:Hz')
ylabel('Power Spectrum, unit:dB')
legend(sprintf('PSD of AR-%d', MMSE_order_y), 'Estimated Spectrum of Noisy Speech'); hold off
