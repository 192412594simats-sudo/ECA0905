
clear; clc; close all;

Fs = 8000;
fc = 1500;
wc = fc / (Fs/2);   % Normalized cutoff frequency
N = 50;             % Filter order

% Window Functions
w_rect  = rectwin(N+1);
w_hamm  = hamming(N+1);
w_hann  = hann(N+1);
w_black = blackman(N+1);

% FIR Filter Design using fir1
b_rect  = fir1(N, wc, w_rect);
b_hamm  = fir1(N, wc, w_hamm);
b_hann  = fir1(N, wc, w_hann);
b_black = fir1(N, wc, w_black);

% Frequency Response Calculation
[H_rect, f] = freqz(b_rect, 1, 1024, Fs);
[H_hamm, ~] = freqz(b_hamm, 1, 1024, Fs);
[H_hann, ~] = freqz(b_hann, 1, 1024, Fs);
[H_black, ~] = freqz(b_black, 1, 1024, Fs);

% Display Numerical Outputs
fprintf('=== EXPERIMENT 4: FIR WINDOWING TECHNIQUES ===\n\n');
fprintf('Filter Order (N): %d\n', N);
fprintf('Cutoff Frequency: %d Hz\n\n', fc);

fprintf('--- RECTANGULAR WINDOW COEFFICIENTS (First 5) ---\n');
disp(b_rect(1:5));

fprintf('--- HAMMING WINDOW COEFFICIENTS (First 5) ---\n');
disp(b_hamm(1:5));

fprintf('--- HANNING WINDOW COEFFICIENTS (First 5) ---\n');
disp(b_hann(1:5));

fprintf('--- BLACKMAN WINDOW COEFFICIENTS (First 5) ---\n');
disp(b_black(1:5));

% Plotting
figure('Name', 'FIR Windowing Techniques Comparison');

% Magnitude Response
subplot(2,1,1);

% Avoid log10(0)
plot(f, 20*log10(max(abs(H_rect), 1e-10)), ...
    'k', 'LineWidth', 1);
hold on;

plot(f, 20*log10(max(abs(H_hamm), 1e-10)), ...
    'b', 'LineWidth', 1.5);

plot(f, 20*log10(max(abs(H_hann), 1e-10)), ...
    'g', 'LineWidth', 1.5);

plot(f, 20*log10(max(abs(H_black), 1e-10)), ...
    'r', 'LineWidth', 1.5);

grid on;
legend('Rectangular', 'Hamming', 'Hann', 'Blackman');
title('FIR Low-Pass Filter Magnitude Response (N = 50)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
ylim([-100 10]);

% Phase Response
subplot(2,1,2);
plot(f, unwrap(angle(H_hamm))*180/pi, ...
    'b', 'LineWidth', 1.5);

grid on;
title('Phase Response (Hamming Window Example)');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
