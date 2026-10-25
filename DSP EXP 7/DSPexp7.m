
clear; clc; close all;

Fs = 10000;
fp = 1500;
fs = 2000;
N = 40;

% Normalized frequency band edges
f_bands = [0 fp fs Fs/2] / (Fs/2);
a_bands = [1 1 0 0];

% 1. Parks-McClellan (Equiripple) Design
b_pm = firpm(N, f_bands, a_bands);

% 2. Frequency-Sampling Design
b_fs = fir2(N, f_bands, a_bands);

% Display numerical outputs
fprintf('=== EXPERIMENT 5: PARKS-MCCLELLAN vs FREQUENCY SAMPLING ===\n\n');
fprintf('Filter Order (N): %d\n', N);
fprintf('Passband Edge: %d Hz, Stopband Edge: %d Hz\n\n', fp, fs);

fprintf('--- PARKS-MCCLELLAN COEFFICIENTS (First 10) ---\n');
disp(b_pm(1:10));

fprintf('--- FREQUENCY-SAMPLING COEFFICIENTS (First 10) ---\n');
disp(b_fs(1:10));

% Calculate frequency responses
[H_pm, f] = freqz(b_pm, 1, 2048, Fs);
[H_fs, ~] = freqz(b_fs, 1, 2048, Fs);

% Avoid log10(0)
mag_pm = 20*log10(max(abs(H_pm), 1e-8));
mag_fs = 20*log10(max(abs(H_fs), 1e-8));

% Plot results
figure('Name', 'Parks-McClellan vs Frequency-Sampling');

subplot(2,1,1);
plot(f, mag_pm, 'b', 'LineWidth', 1.5);
hold on;
plot(f, mag_fs, 'r--', 'LineWidth', 1.5);
grid on;

legend('Parks-McClellan (Equiripple)', ...
    'Frequency-Sampling', 'Location', 'best');
title(['FIR Filter Magnitude Response Comparison (Order N = ', ...
    num2str(N), ')']);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([0 Fs/2]);
ylim([-80 10]);

subplot(2,1,2);
plot(f, angle(H_pm)*180/pi, 'b', 'LineWidth', 1.5);
hold on;
plot(f, angle(H_fs)*180/pi, 'r--', 'LineWidth', 1.5);
grid on;

legend('Parks-McClellan', 'Frequency-Sampling', ...
    'Location', 'best');
title('Phase Response Comparison');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
xlim([0 Fs/2]);
ylim([-200 200]);
