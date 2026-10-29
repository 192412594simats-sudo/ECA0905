clc;
clear;
close all;

%% 1. SYNTHETIC SPEECH GENERATION
Fs = 16000;                 % Sampling rate (Hz)
duration = 0.5;             % Duration (seconds)
t = 0:1/Fs:duration-1/Fs;

f1 = 250;
f2 = 700;
x_speech = sin(2*pi*f1*t) + 0.5*sin(2*pi*f2*t);

%% 2. FRAME-BY-FRAME LPC VOCODER SYNTHESIS
p = 12;                     % LPC prediction order
N_frame = 320;
hop = 160;
win = hamming(N_frame).';

num_frames = floor((length(x_speech)-N_frame)/hop) + 1;
x_synth = zeros(1, length(x_speech));
win_gain = zeros(1, length(x_speech));

rng(42);                    % Fixed random seed

for m = 1:num_frames
    idx = (m-1)*hop + (1:N_frame);
    frame = x_speech(idx) .* win;

    % LPC coefficients and prediction-error gain
    [a, g] = lpc(frame, p);

    % White-noise excitation and LPC synthesis filter
    excitation = sqrt(max(g, 0)) * randn(1, N_frame);
    synth_frame = filter(1, a, excitation);

    % Overlap-add with window normalization
    x_synth(idx) = x_synth(idx) + synth_frame .* win;
    win_gain(idx) = win_gain(idx) + win.^2;
end

valid_idx = win_gain > 1e-8;
x_synth(valid_idx) = x_synth(valid_idx) ./ win_gain(valid_idx);

%% 3. SINGLE-FRAME SPECTRAL ANALYSIS
frame_start = 161;
frame_idx = frame_start:(frame_start + N_frame - 1);
frame_sample = x_speech(frame_idx) .* win;

[a_sample, g_sample] = lpc(frame_sample, p);

[H_lpc, f_lpc] = freqz(sqrt(max(g_sample, 0)), a_sample, 1024, Fs);
[Pxx, f_pwelch] = pwelch(frame_sample, win, 0, 1024, Fs);

lpc_db = 20*log10(max(abs(H_lpc), 1e-12));
pwelch_db = 10*log10(max(Pxx, 1e-12));

%% 4. PLOTS AND EXPORT
figW = 8.5;
figH = 5.5;

% Figure 20A: Original signal
figA = figure('Units', 'centimeters', ...
    'Position', [2, 2, figW, figH]);

plot(t(1:640), x_speech(1:640), 'b', 'LineWidth', 1.2);
grid on;
box on;
title('Original Synthetic Speech Signal');
xlabel('Time (s)');
ylabel('Amplitude');

exportgraphics(figA, 'Exp20_SubfigA_Original.png', 'Resolution', 600);

% Figure 20B: LPC synthesized signal
figB = figure('Units', 'centimeters', ...
    'Position', [11, 2, figW, figH]);

plot(t(1:640), x_synth(1:640), 'r', 'LineWidth', 1.2);
grid on;
box on;
title('LPC Vocoder Synthesized Signal');
xlabel('Time (s)');
ylabel('Amplitude');

exportgraphics(figB, 'Exp20_SubfigB_Synthesized.png', 'Resolution', 600);

% Figure 20C: LPC envelope and signal spectrum
figC = figure('Units', 'centimeters', ...
    'Position', [20, 2, figW, figH]);

plot(f_pwelch, pwelch_db, 'b', 'LineWidth', 1.2);
hold on;
plot(f_lpc, lpc_db, 'r', 'LineWidth', 1.2);
grid on;
box on;
xlim([0 Fs/2]);
title('Single-Frame Spectrum and LPC Envelope');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Welch Spectrum', 'LPC Envelope');

exportgraphics(figC, 'Exp20_SubfigC_Spectrum.png', 'Resolution', 600);