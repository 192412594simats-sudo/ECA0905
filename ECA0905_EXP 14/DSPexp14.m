clear;
clc;
close all;

%% 1. PARAMETER DEFINITIONS
fc_carrier = 28e9;       % Carrier frequency (28 GHz; informational only)
Fs_in = 1000e6;          % Input sampling frequency (1000 MHz)
M = 4;                   % Decimation factor
L = 2;                   % Interpolation factor
Fs_dec = Fs_in / M;      % Decimated sampling frequency (250 MHz)
Fs_out = Fs_dec * L;     % Output sampling frequency (500 MHz)
N_samples = 4000;        % Number of input samples

t_in = (0:N_samples-1) / Fs_in;

% Generate baseband multitone signal
f1 = 15e6;
f2 = 35e6;
x_in = sin(2*pi*f1*t_in) + 0.5*cos(2*pi*f2*t_in);

%% 2. DECIMATION: ANTI-ALIASING FILTER AND DOWNSAMPLING
N_fir = 60;

% Cutoff at the new Nyquist frequency, normalized to input Nyquist
fc_dec = (Fs_dec/2) / (Fs_in/2);
b_dec = fir1(N_fir, fc_dec, 'low');

x_filtered_dec = filter(b_dec, 1, x_in);
x_dec = x_filtered_dec(1:M:end);

%% 3. INTERPOLATION: UPSAMPLING AND ANTI-IMAGING FILTER
x_up = zeros(1, length(x_dec)*L);
x_up(1:L:end) = x_dec;

% Cutoff at the decimated Nyquist frequency, normalized to output Nyquist
fc_interp = (Fs_dec/2) / (Fs_out/2);
b_interp = fir1(N_fir, fc_interp, 'low');

x_out = L * filter(b_interp, 1, x_up);

%% 4. DECIMATION FILTER PERFORMANCE
[H, f_response] = freqz(b_dec, 1, 4096, Fs_in);
H_mag_db = 20*log10(max(abs(H), 1e-12));

pb_idx = f_response <= 0.8*(Fs_dec/2);
sb_idx = f_response >= Fs_dec/2;

passband_ripple = max(H_mag_db(pb_idx)) - min(H_mag_db(pb_idx));
stopband_atten = -max(H_mag_db(sb_idx));

%% 5. COMMAND WINDOW RESULTS
fprintf('========================================================\n');
fprintf(' EXPERIMENT 14: MULTIRATE BASEBAND PROCESSING RESULTS\n');
fprintf('========================================================\n');
fprintf('Carrier frequency (reference) : %.2f GHz\n', fc_carrier/1e9);
fprintf('Input sampling rate            : %.2f MHz\n', Fs_in/1e6);
fprintf('Decimation factor (M)           : %d\n', M);
fprintf('Downsampled sampling rate       : %.2f MHz\n', Fs_dec/1e6);
fprintf('Interpolation factor (L)        : %d\n', L);
fprintf('Output sampling rate            : %.2f MHz\n\n', Fs_out/1e6);

fprintf('Sequence lengths:\n');
fprintf(' Input                         : %d samples\n', length(x_in));
fprintf(' Decimated                     : %d samples\n', length(x_dec));
fprintf(' Interpolated output           : %d samples\n\n', length(x_out));

fprintf('Decimation filter passband ripple: %.2f dB\n', passband_ripple);
fprintf('Decimation filter stopband attenuation: %.2f dB\n\n', stopband_atten);
fprintf('Status: Multirate filtering completed.\n');
fprintf('Note: The generated signal is baseband; no 28 GHz carrier is applied.\n');
fprintf('========================================================\n');

%% 6. SIGNAL SPECTRA
Nfft = 8192;

[X_in, f_in] = one_sided_spectrum(x_in, Fs_in, Nfft);
[X_dec, f_dec] = one_sided_spectrum(x_dec, Fs_dec, Nfft);
[X_out, f_out] = one_sided_spectrum(x_out, Fs_out, Nfft);

figure('Name', 'Experiment 14: Multirate Processing', ...
       'NumberTitle', 'off');

subplot(3,1,1);
plot(f_in/1e6, X_in, 'b', 'LineWidth', 1.2);
grid on;
title('Original Baseband Signal Spectrum (Fs = 1000 MHz)');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
xlim([0 Fs_in/2/1e6]);
ylim([-80 5]);

subplot(3,1,2);
plot(f_dec/1e6, X_dec, 'r', 'LineWidth', 1.2);
grid on;
title('Decimated Signal Spectrum (M = 4, Fs = 250 MHz)');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
xlim([0 Fs_dec/2/1e6]);
ylim([-80 5]);

subplot(3,1,3);
plot(f_out/1e6, X_out, 'm', 'LineWidth', 1.2);
grid on;
title('Interpolated Output Spectrum (L = 2, Fs = 500 MHz)');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
xlim([0 Fs_out/2/1e6]);
ylim([-80 5]);

%% LOCAL FUNCTION: ONE-SIDED, NORMALIZED SPECTRUM
function [mag_db, f] = one_sided_spectrum(x, Fs, Nfft)
    X = fft(x, Nfft);
    nPositive = floor(Nfft/2) + 1;

    mag = abs(X(1:nPositive));
    mag = mag / max(max(mag), eps);
    mag_db = 20*log10(max(mag, 1e-12));

    f = (0:nPositive-1) * (Fs/Nfft);
end