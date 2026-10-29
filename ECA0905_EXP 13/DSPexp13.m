clc;
clear;
close all;

%% 1. INPUT PARAMETERS
N = 24;                     % FIR filter order
M = 2;                      % Decimation factor
Fs = 8000;                  % Sampling frequency (Hz)
N_samples = 128;            % Number of input samples

% FIR prototype filter
b = fir1(N, 0.4);

% Polyphase decomposition
E0 = b(1:M:end);            % Even-indexed coefficients
E1 = b(2:M:end);            % Odd-indexed coefficients

% Input signal: 500 Hz sinusoid
t = (0:N_samples-1) / Fs;
x = sin(2*pi*500*t);

%% 2. CONVENTIONAL AND POLYPHASE DECIMATION

% Direct path: filter, then downsample
y_direct_full = filter(b, 1, x);
y_direct = y_direct_full(1:M:end);

% Polyphase path: split input phases and filter at low rate
x0 = x(1:M:end);
x_delayed = [0, x(1:end-1)];
x1 = x_delayed(1:M:end);

y_poly = filter(E0, 1, x0) + filter(E1, 1, x1);

% Verification
rec_error = max(abs(y_direct - y_poly));

%% 3. COMPUTATIONAL COMPLEXITY
N_taps = length(b);
N_out = length(y_direct);

% Count one multiplication per tap per output sample
MAC_direct = N_out * N_taps;
MAC_poly = N_out * (length(E0) + length(E1));
comp_saving = (1 - MAC_poly / MAC_direct) * 100;

%% 4. COMMAND WINDOW DISPLAY
fprintf('=============================================================\n');
fprintf(' EXPERIMENT 13: POLYPHASE DECOMPOSITION OF DECIMATION FILTER\n');
fprintf('=============================================================\n');
fprintf('Decimation factor (M)       : %d\n', M);
fprintf('Filter order (N)            : %d (%d taps)\n', N, N_taps);
fprintf('Input signal                : 500 Hz sinusoid (%d samples)\n', N_samples);
fprintf('Sampling frequency (Fs)     : %d Hz\n\n', Fs);

fprintf('Polyphase sub-filter coefficients:\n');
fprintf(' E0(z) = [');
fprintf('%.4f ', E0);
fprintf(']\n');

fprintf(' E1(z) = [');
fprintf('%.4f ', E1);
fprintf(']\n\n');

fprintf('COMPUTATIONAL COMPLEXITY:\n');
fprintf(' Direct implementation     : %d MACs\n', MAC_direct);
fprintf(' Polyphase implementation  : %d MACs\n', MAC_poly);
fprintf(' Computational savings     : %.2f%%\n\n', comp_saving);

fprintf('VERIFICATION:\n');
fprintf(' Maximum absolute difference: %.4e\n', rec_error);

if rec_error < 1e-12
    fprintf(' Status: Polyphase output matches direct output.\n');
else
    fprintf(' Status: Outputs differ; check phase alignment.\n');
end

fprintf('=============================================================\n');

%% 5. FIGURES
figW = 8.5;
figH = 5.5;

% Figure 13A: Prototype impulse response
figA = figure('Units', 'centimeters', ...
    'Position', [2, 2, figW, figH]);

stem(0:N, b, 'filled', 'LineWidth', 1.2);
grid on;
box on;
title('Prototype FIR Filter Impulse Response');
xlabel('Sample Index (n)');
ylabel('Coefficient Value');

exportgraphics(figA, 'Exp13_SubfigA_Prototype.png', 'Resolution', 600);

% Figure 13B: Polyphase sub-filter coefficients
figB = figure('Units', 'centimeters', ...
    'Position', [11, 2, figW, figH]);

stem(0:length(E0)-1, E0, 'b', 'filled', 'LineWidth', 1.2);
hold on;
stem(0:length(E1)-1, E1, 'r--', 'LineWidth', 1.2);
grid on;
box on;
title('Polyphase Sub-Filter Coefficients (M = 2)');
xlabel('Sub-Filter Sample Index');
ylabel('Coefficient Value');
legend('Branch E0[n]', 'Branch E1[n]', 'FontSize', 7);

exportgraphics(figB, 'Exp13_SubfigB_Polyphase.png', 'Resolution', 600);

% Figure 13C: Compare outputs
figC = figure('Units', 'centimeters', ...
    'Position', [20, 2, figW, figH]);

numPlot = min(50, length(y_direct));
plot(0:numPlot-1, y_direct(1:numPlot), 'b', 'LineWidth', 1.5);
hold on;
plot(0:numPlot-1, y_poly(1:numPlot), 'r--', 'LineWidth', 1.5);
grid on;
box on;
title('Time-Domain Output Comparison');
xlabel('Output Sample Index');
ylabel('Amplitude');
legend('Direct Decimation', 'Polyphase Decimation', 'FontSize', 7);

exportgraphics(figC, 'Exp13_SubfigC_OutputComparison.png', 'Resolution', 600);