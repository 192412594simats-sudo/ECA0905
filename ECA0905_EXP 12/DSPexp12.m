clc;
clear;
close all;

%% 1. SIGNAL GENERATION
Fs = 8000;                  % Original sampling frequency (Hz)
f1 = 500;
f2 = 1500;
M = 2;                      % Decimation factor
L = 2;                      % Interpolation factor
N = 1000;                    % Number of input samples

fc_dec = 1800;              % Decimation filter cutoff (Hz)
fc_interp = 3600;           % Interpolation filter cutoff (Hz)

n = 0:N-1;
t = n / Fs;
x = sin(2*pi*f1*t) + 0.8*sin(2*pi*f2*t);

%% 2. DECIMATION
Fs_dec = Fs / M;

% Direct downsampling (no anti-aliasing filter)
x_down_direct = x(1:M:end);

% Anti-aliasing filter, designed at the original sampling frequency
b_dec = fir1(60, fc_dec/(Fs/2));
x_dec_filtered = filter(b_dec, 1, x);
x_dec = x_dec_filtered(1:M:end);

%% 3. INTERPOLATION
Fs_interp = Fs * L;

% Insert L-1 zeros between input samples
x_up = zeros(1, length(x)*L);
x_up(1:L:end) = x;

% Anti-imaging filter, designed at the interpolated sampling frequency
b_interp = L * fir1(60, fc_interp/(Fs_interp/2));
x_interp_filtered = filter(b_interp, 1, x_up);

%% 4. DISPLAY RESULTS
fprintf('\n==================================================\n');
fprintf(' EXPERIMENT 12: MULTIRATE SIGNAL PROCESSING\n');
fprintf('==================================================\n');
fprintf('Original sampling frequency       : %d Hz\n', Fs);
fprintf('Decimation factor                 : %d\n', M);
fprintf('Sampling frequency after decimation: %d Hz\n', Fs_dec);
fprintf('Interpolation factor              : %d\n', L);
fprintf('Sampling frequency after interpolation: %d Hz\n', Fs_interp);
fprintf('\nSignal lengths:\n');
fprintf('Original signal                   : %d samples\n', length(x));
fprintf('Decimated signal                  : %d samples\n', length(x_dec));
fprintf('Interpolated signal               : %d samples\n', length(x_interp_filtered));
fprintf('\nStatus: Decimation and interpolation completed.\n');
fprintf('==================================================\n');

%% 5. SPECTRUM COMPUTATION
Nfft = 2048;

[X_mag, f_x] = compute_spectrum(x, Fs, Nfft);
[X_down_mag, f_down] = compute_spectrum(x_down_direct, Fs_dec, Nfft);
[X_dec_mag, f_dec] = compute_spectrum(x_dec, Fs_dec, Nfft);
[X_up_mag, f_up] = compute_spectrum(x_up, Fs_interp, Nfft);
[X_interp_mag, f_interp] = compute_spectrum(x_interp_filtered, Fs_interp, Nfft);

%% 6. ORIGINAL SIGNAL PLOTS
figA = figure('Units', 'centimeters', 'Position', [2, 2, 8.5, 6.5]);

subplot(2,1,1);
idx = 1:min(100, length(x));
plot(t(idx)*1000, x(idx), 'b', 'LineWidth', 1.2);
grid on;
title('Original Input Waveform');
xlabel('Time (ms)');
ylabel('Amplitude');

subplot(2,1,2);
plot(f_x, X_mag, 'b', 'LineWidth', 1.2);
grid on;
title('Original Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

exportgraphics(figA, 'Exp12_SubfigA_Original.png', 'Resolution', 600);

%% 7. DECIMATED SIGNAL PLOTS
figB = figure('Units', 'centimeters', 'Position', [11, 2, 8.5, 6.5]);

t_dec = (0:length(x_dec)-1) / Fs_dec;
idx = 1:min(50, length(x_dec));

subplot(2,1,1);
plot(t_dec(idx)*1000, x_dec(idx), 'r', 'LineWidth', 1.2);
grid on;
title('Filtered Decimated Waveform (M = 2)');
xlabel('Time (ms)');
ylabel('Amplitude');

subplot(2,1,2);
plot(f_down, X_down_mag, 'k--', 'LineWidth', 1);
hold on;
plot(f_dec, X_dec_mag, 'r', 'LineWidth', 1.2);
grid on;
title('Decimated Signal Spectra');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Direct Downsampling', 'Filtered Decimation');

exportgraphics(figB, 'Exp12_SubfigB_Decimated.png', 'Resolution', 600);

%% 8. INTERPOLATED SIGNAL PLOTS
figC = figure('Units', 'centimeters', 'Position', [20, 2, 8.5, 6.5]);

t_interp = (0:length(x_interp_filtered)-1) / Fs_interp;
idx = 1:min(200, length(x_interp_filtered));

subplot(2,1,1);
plot(t_interp(idx)*1000, x_interp_filtered(idx), 'm', 'LineWidth', 1.2);
grid on;
title('Interpolated Waveform (L = 2)');
xlabel('Time (ms)');
ylabel('Amplitude');

subplot(2,1,2);
plot(f_up, X_up_mag, 'k--', 'LineWidth', 1);
hold on;
plot(f_interp, X_interp_mag, 'm', 'LineWidth', 1.2);
grid on;
title('Interpolated Signal Spectra');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Zero Insertion', 'Filtered Interpolation');

exportgraphics(figC, 'Exp12_SubfigC_Interpolated.png', 'Resolution', 600);

%% HELPER FUNCTION
function [X_mag, f] = compute_spectrum(signal, Fs_curr, Nfft)
    X = fft(signal, Nfft);
    X_mag = abs(fftshift(X)) / length(signal);
    f = (-Nfft/2 : Nfft/2-1) * (Fs_curr/Nfft);
end